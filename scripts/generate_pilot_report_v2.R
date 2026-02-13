#!/usr/bin/env Rscript
# =============================================================================
# generate_pilot_report_v2.R — Generate calibration summary + pilot figures
# for IRW-based pilot (v2)
# =============================================================================

library(nimble)
library(DPMirt)
library(IRTsimrel)
library(tidyverse)
library(yaml)
library(digest)

source("R/00_config.R")
source("R/01_calibrate.R")
source("R/05_postprocess.R")

cfg <- load_config()

# =============================================================================
# 1. Rebuild calibration summary from cached .rds files
# =============================================================================
cat("=== Rebuilding calibration summary from IRW calibration ===\n")
calib_results <- run_all_calibrations(cfg)  # will use cache
summary_tbl <- build_calibration_summary(calib_results)

# Add latent_shape column by parsing the cell ID
summary_tbl$latent_shape <- sapply(strsplit(summary_tbl$calib_cell_id, "_"), function(x) {
  # Format: model_shape_rhoNN
  # rasch_normal_rho50 -> normal
  # 2pl_bimodal_rho60 -> bimodal
  # rasch_skew_pos_rho70 -> skew_pos
  parts <- x
  # Find where "rho" starts
  rho_idx <- grep("^rho", parts)
  model_part <- parts[1]
  shape_parts <- parts[2:(rho_idx - 1)]
  paste(shape_parts, collapse = "_")
})

saveRDS(summary_tbl, file.path(cfg$paths$calibration_dir, "calibration_summary.rds"))
write.csv(summary_tbl, file.path(cfg$paths$calibration_dir, "calibration_summary.csv"),
          row.names = FALSE)
cat(sprintf("Calibration summary: %d cells saved\n", nrow(summary_tbl)))
cat(sprintf("  Max rho error: %.6e\n", max(summary_tbl$rho_error)))

# =============================================================================
# 2. Calibration figure (IRW)
# =============================================================================
cat("\n=== Generating calibration figure (IRW) ===\n")

calib_plot_df <- summary_tbl %>%
  mutate(
    model_label = ifelse(irt_model == "rasch", "Rasch", "2PL"),
    shape_label = dplyr::case_when(
      latent_shape == "normal"   ~ "Normal",
      latent_shape == "bimodal"  ~ "Bimodal",
      latent_shape == "skew_pos" ~ "Skewed",
      TRUE ~ latent_shape
    )
  )

p_calib <- ggplot(calib_plot_df,
                  aes(x = target_rho, y = c_star, color = shape_label, shape = shape_label)) +
  geom_point(size = 3) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ model_label) +
  scale_color_manual(values = c("Normal" = "#2166AC", "Bimodal" = "#B2182B", "Skewed" = "#4DAF4A")) +
  labs(
    title = "EQC Calibration: Scaling Constant c* vs Target Reliability (IRW Items)",
    subtitle = sprintf("All %d cells converged | Max rho error: %.2e | Item source: IRW",
                        nrow(summary_tbl), max(summary_tbl$rho_error)),
    x = expression("Target Reliability " * rho),
    y = expression("Scaling Constant " * c * "*"),
    color = "Latent Shape",
    shape = "Latent Shape"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

ggsave(file.path(cfg$paths$figures_dir, "phase1a_calibration_c_star_irw.png"),
       p_calib, width = 10, height = 5, dpi = 150)
cat("Saved: phase1a_calibration_c_star_irw.png\n")

# =============================================================================
# 3. Pilot results figures (IRW + parallel)
# =============================================================================
cat("\n=== Generating pilot figures (IRW) ===\n")

summary_df <- readRDS(file.path(cfg$paths$results_dir, "pilot_results_summary.rds"))

# RMSE by condition
p_rmse <- plot_pilot_rmse(
  summary_df,
  save_path = file.path(cfg$paths$figures_dir, "pilot_rmse_by_condition_irw.png")
)
cat("Saved: pilot_rmse_by_condition_irw.png\n")

# KS by condition
p_ks <- plot_pilot_ks(
  summary_df,
  save_path = file.path(cfg$paths$figures_dir, "pilot_ks_by_condition_irw.png")
)
cat("Saved: pilot_ks_by_condition_irw.png\n")

# Overall RMSE by method (pooled across conditions)
overall_rmse <- summary_df %>%
  filter(!is.na(mean_rmse)) %>%
  group_by(prior, method) %>%
  summarise(overall_rmse = mean(mean_rmse), .groups = "drop") %>%
  mutate(
    prior_label = dplyr::case_when(
      prior == "gaussian"   ~ "Gaussian",
      prior == "dp_broad"   ~ "DP-Broad",
      prior == "dp_focused" ~ "DP-Focused",
      TRUE ~ prior
    ),
    method_label = toupper(method)
  )

p_overall_rmse <- ggplot(overall_rmse, aes(x = method_label, y = overall_rmse, fill = prior_label)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("Gaussian" = "#2166AC", "DP-Broad" = "#B2182B", "DP-Focused" = "#4DAF4A")) +
  labs(
    title = "Pilot Overall RMSE by Method and Prior (IRW Items)",
    x = "Summary Method", y = "Mean RMSE (pooled across 4 conditions)", fill = "Prior"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

ggsave(file.path(cfg$paths$figures_dir, "pilot_overall_rmse_irw.png"),
       p_overall_rmse, width = 8, height = 5, dpi = 150)
cat("Saved: pilot_overall_rmse_irw.png\n")

# Overall KS by method
overall_ks <- summary_df %>%
  filter(!is.na(mean_ks)) %>%
  group_by(prior, method) %>%
  summarise(overall_ks = mean(mean_ks), .groups = "drop") %>%
  mutate(
    prior_label = dplyr::case_when(
      prior == "gaussian"   ~ "Gaussian",
      prior == "dp_broad"   ~ "DP-Broad",
      prior == "dp_focused" ~ "DP-Focused",
      TRUE ~ prior
    ),
    method_label = toupper(method)
  )

p_overall_ks <- ggplot(overall_ks, aes(x = method_label, y = overall_ks, fill = prior_label)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("Gaussian" = "#2166AC", "DP-Broad" = "#B2182B", "DP-Focused" = "#4DAF4A")) +
  labs(
    title = "Pilot Overall KS Distance by Method and Prior (IRW Items)",
    x = "Summary Method", y = "Mean KS Distance (pooled across 4 conditions)", fill = "Prior"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

ggsave(file.path(cfg$paths$figures_dir, "pilot_overall_ks_irw.png"),
       p_overall_ks, width = 8, height = 5, dpi = 150)
cat("Saved: pilot_overall_ks_irw.png\n")

# =============================================================================
# 4. Comparison: parametric vs IRW item parameters (illustrative)
# =============================================================================
cat("\n=== Generating IRW item parameter distribution figure ===\n")

# Collect beta and lambda from all 30 calibration cells
item_df <- bind_rows(lapply(names(calib_results), function(cid) {
  res <- calib_results[[cid]]
  data.frame(
    cell_id = cid,
    model = res$model,
    item = seq_along(res$beta_vec),
    beta = res$beta_vec,
    lambda_scaled = res$lambda_scaled,
    stringsAsFactors = FALSE
  )
}))

p_items <- ggplot(item_df, aes(x = beta)) +
  geom_histogram(bins = 40, fill = "#2166AC", alpha = 0.7, color = "white") +
  facet_wrap(~ ifelse(model == "rasch", "Rasch", "2PL"), scales = "free_y") +
  labs(
    title = "IRW Item Difficulty Distribution (Across All 30 Calibration Cells)",
    subtitle = sprintf("Total items: %d | Source: Item Response Warehouse", nrow(item_df)),
    x = expression("Item Difficulty " * beta),
    y = "Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(cfg$paths$figures_dir, "irw_item_difficulty_distribution.png"),
       p_items, width = 10, height = 4.5, dpi = 150)
cat("Saved: irw_item_difficulty_distribution.png\n")

cat("\n=== All figures generated ===\n")
