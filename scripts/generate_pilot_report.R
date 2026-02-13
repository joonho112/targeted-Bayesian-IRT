#!/usr/bin/env Rscript
# =============================================================================
# generate_pilot_report.R  —  Generate figures and report from pilot results
# =============================================================================
# Run AFTER test_pilot.R completes. Reads saved results and generates:
#   1. RMSE bar chart by condition x prior x method
#   2. KS bar chart by condition x prior x method
#   3. MSELR bar chart
#   4. Theoretical validation table (PM best MSE, GR best KS)
#   5. Timing projection table
#   6. Markdown progress report
# =============================================================================

library(tidyverse)
library(yaml)

source("R/00_config.R")
source("R/05_postprocess.R")

cfg <- load_config()

# =============================================================================
# Load pilot results
# =============================================================================
raw_path <- file.path(cfg$paths$results_dir, "pilot_results_raw.rds")
sum_path <- file.path(cfg$paths$results_dir, "pilot_results_summary.rds")

if (!file.exists(raw_path)) {
  stop("Pilot results not found. Run test_pilot.R first.")
}

results_df <- readRDS(raw_path)
summary_df <- readRDS(sum_path)

# Focus on theta
theta_df <- results_df %>% filter(parameter == "theta")

cat(sprintf("Loaded %d result rows (%d theta)\n", nrow(results_df), nrow(theta_df)))

# =============================================================================
# Figure 1: RMSE by condition
# =============================================================================
fig1_path <- file.path(cfg$paths$figures_dir, "pilot_rmse_by_condition.png")
plot_pilot_rmse(summary_df, save_path = fig1_path)
cat(sprintf("Saved: %s\n", fig1_path))

# =============================================================================
# Figure 2: KS by condition
# =============================================================================
fig2_path <- file.path(cfg$paths$figures_dir, "pilot_ks_by_condition.png")
plot_pilot_ks(summary_df, save_path = fig2_path)
cat(sprintf("Saved: %s\n", fig2_path))

# =============================================================================
# Figure 3: Method comparison (PM vs CB vs GR) pooled across conditions
# =============================================================================
method_summary <- theta_df %>%
  mutate(rmse = sqrt(msel)) %>%
  group_by(prior, method) %>%
  summarise(
    mean_rmse  = mean(rmse, na.rm = TRUE),
    mean_mselr = mean(mselr, na.rm = TRUE),
    mean_ks    = mean(ks, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    prior_label = case_when(
      prior == "gaussian" ~ "Gaussian",
      prior == "dp_broad" ~ "DP-Broad",
      prior == "dp_focused" ~ "DP-Focused"
    ),
    method_label = toupper(method)
  )

# RMSE comparison
p3 <- ggplot(method_summary, aes(x = method_label, y = mean_rmse, fill = prior_label)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("Gaussian"="#2166AC", "DP-Broad"="#B2182B", "DP-Focused"="#4DAF4A")) +
  labs(title = "Overall Mean RMSE by Prior and Method",
       subtitle = "Averaged over all pilot conditions and replications",
       x = "Summary Method", y = "Mean RMSE", fill = "Prior") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

fig3_path <- file.path(cfg$paths$figures_dir, "pilot_overall_rmse.png")
ggsave(fig3_path, p3, width = 8, height = 5, dpi = 150)
cat(sprintf("Saved: %s\n", fig3_path))

# KS comparison
p4 <- ggplot(method_summary, aes(x = method_label, y = mean_ks, fill = prior_label)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("Gaussian"="#2166AC", "DP-Broad"="#B2182B", "DP-Focused"="#4DAF4A")) +
  labs(title = "Overall Mean KS Distance by Prior and Method",
       subtitle = "Theory predicts GR has lowest KS (best distribution recovery)",
       x = "Summary Method", y = "Mean KS Distance", fill = "Prior") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

fig4_path <- file.path(cfg$paths$figures_dir, "pilot_overall_ks.png")
ggsave(fig4_path, p4, width = 8, height = 5, dpi = 150)
cat(sprintf("Saved: %s\n", fig4_path))

# =============================================================================
# Sanity check table
# =============================================================================
cat("\n=== THEORETICAL VALIDATION ===\n\n")

# For each condition x prior, check PM has best RMSE and GR has best KS
validation <- theta_df %>%
  mutate(rmse = sqrt(msel)) %>%
  group_by(irt_model, latent_shape, N, target_rho, prior, method) %>%
  summarise(mean_rmse = mean(rmse), mean_ks = mean(ks), .groups = "drop") %>%
  group_by(irt_model, latent_shape, N, target_rho, prior) %>%
  summarise(
    best_rmse_method = method[which.min(mean_rmse)],
    best_ks_method   = method[which.min(mean_ks)],
    .groups = "drop"
  )

cat("PM should have best (lowest) RMSE:\n")
pm_pct <- mean(validation$best_rmse_method == "pm") * 100
cat(sprintf("  PM is best in %.0f%% of condition-prior combos (%d/%d)\n",
    pm_pct, sum(validation$best_rmse_method == "pm"), nrow(validation)))

cat("\nGR should have best (lowest) KS:\n")
gr_pct <- mean(validation$best_ks_method == "gr") * 100
cat(sprintf("  GR is best in %.0f%% of condition-prior combos (%d/%d)\n",
    gr_pct, sum(validation$best_ks_method == "gr"), nrow(validation)))

# =============================================================================
# Timing table
# =============================================================================
cat("\n=== TIMING PROJECTION ===\n\n")

timing <- theta_df %>%
  distinct(condition_id, irt_model, N, rep_id, prior, .keep_all = TRUE) %>%
  group_by(irt_model, N) %>%
  summarise(
    mean_time = mean(time_secs, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  )

cat(sprintf("%-8s  %5s  %8s\n", "Model", "N", "sec/fit"))
cat(paste(rep("-", 25), collapse = ""), "\n")
for (i in 1:nrow(timing)) {
  cat(sprintf("%-8s  %5d  %8.1f\n", timing$irt_model[i], timing$N[i], timing$mean_time[i]))
}

# Full simulation projection
mean_time_per_fit <- mean(theta_df$time_secs, na.rm = TRUE) / 3  # 3 methods per fit
total_fits <- 120 * 100 * 3  # conditions * reps * priors
total_hours_seq <- (mean_time_per_fit * total_fits) / 3600
cat(sprintf("\nEstimated full simulation (sequential): %.0f hours\n", total_hours_seq))
cat(sprintf("With 4 parallel workers: ~%.0f hours\n", total_hours_seq / 4))
cat(sprintf("With 8 parallel workers: ~%.0f hours\n", total_hours_seq / 8))

cat("\nReport generation complete.\n")
