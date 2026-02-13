#!/usr/bin/env Rscript
# =============================================================================
# test_phase1a_calibration.R
# =============================================================================
# Run EQC calibration for all 30 cells and verify results.
# Generates: calibration summary table + diagnostic figure
# =============================================================================

library(IRTsimrel)
library(tidyverse)
library(yaml)
library(digest)

source("R/00_config.R")
source("R/01_calibrate.R")

cfg <- load_config()

# Clear old cache for fresh run
old_files <- list.files(cfg$paths$calibration_dir, pattern = "^eqc_.*\\.rds$", full.names = TRUE)
if (length(old_files) > 0) {
  file.remove(old_files)
  cat(sprintf("Cleared %d old calibration files.\n", length(old_files)))
}

# =============================================================================
# Run all 30 calibrations
# =============================================================================
cat("\n====================================================\n")
cat("  Phase 1a: EQC Calibration (30 cells)\n")
cat("====================================================\n\n")

t0 <- Sys.time()
calib_results <- run_all_calibrations(cfg)
t1 <- Sys.time()
elapsed <- as.numeric(difftime(t1, t0, units = "secs"))
cat(sprintf("\nTotal calibration time: %.1f seconds\n\n", elapsed))

# =============================================================================
# Build summary and verify
# =============================================================================
summary_df <- build_calibration_summary(calib_results)

cat("====================================================\n")
cat("  CALIBRATION RESULTS (30 cells)\n")
cat("====================================================\n\n")

# Print formatted table
cat(sprintf("%-28s  %7s  %6s  %8s  %10s  %8s\n",
    "Cell ID", "c*", "target", "achieved", "error", "status"))
cat(paste(rep("-", 80), collapse = ""), "\n")

for (i in 1:nrow(summary_df)) {
  r <- summary_df[i, ]
  status_mark <- ifelse(r$root_status == "uniroot_success", "OK", "FAIL")
  cat(sprintf("%-28s  %7.4f  %6.2f  %8.6f  %10.8f  %8s\n",
      r$calib_cell_id, r$c_star, r$target_rho, r$achieved_rho,
      r$rho_error, status_mark))
}

# =============================================================================
# Verification checks
# =============================================================================
cat("\n====================================================\n")
cat("  VERIFICATION\n")
cat("====================================================\n\n")

all_converged <- all(summary_df$root_status == "uniroot_success")
all_within_tol <- all(summary_df$rho_error < 1e-4)
max_error <- max(summary_df$rho_error)

cat(sprintf("Total cells:      %d / 30 expected\n", nrow(summary_df)))
cat(sprintf("All converged:    %s\n", ifelse(all_converged, "PASS", "FAIL")))
cat(sprintf("Max rho error:    %.8f\n", max_error))
cat(sprintf("Within tol:       %s (tol = 1e-4)\n", ifelse(all_within_tol, "PASS", "FAIL")))
cat(sprintf("c* range:         [%.4f, %.4f]\n", min(summary_df$c_star), max(summary_df$c_star)))
cat(sprintf("c* within [0.3, 3.0]: %s\n",
    ifelse(all(summary_df$c_star >= 0.3 & summary_df$c_star <= 3.0), "PASS", "WARN")))

# =============================================================================
# Generate diagnostic figure: c* vs target reliability by model and shape
# =============================================================================
cat("\nGenerating calibration diagnostic figure...\n")

# Parse cell_id into components
summary_df <- summary_df %>%
  mutate(
    latent_shape = case_when(
      grepl("_normal_", calib_cell_id) ~ "Normal",
      grepl("_bimodal_", calib_cell_id) ~ "Bimodal",
      grepl("_skew_pos_", calib_cell_id) ~ "Skewed"
    )
  )

p <- ggplot(summary_df, aes(x = target_rho, y = c_star,
                             color = latent_shape, shape = latent_shape)) +
  geom_point(size = 3) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ irt_model, labeller = labeller(irt_model = c("rasch" = "Rasch", "2pl" = "2PL"))) +
  scale_x_continuous(breaks = seq(0.5, 0.9, 0.1)) +
  scale_color_manual(values = c("Normal" = "#2166AC", "Bimodal" = "#B2182B", "Skewed" = "#4DAF4A")) +
  labs(
    title = "EQC Calibration Results: Scaling Factor c* by Target Reliability",
    subtitle = sprintf("I = 25 items, M = %d quadrature points, %d cells total",
                       cfg$calibration$M_quadrature, nrow(summary_df)),
    x = expression("Target Reliability " * bar(omega)),
    y = expression("Discrimination Scale " * c * "*"),
    color = "Latent Shape",
    shape = "Latent Shape"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

fig_path <- file.path(cfg$paths$figures_dir, "phase1a_calibration_c_star.png")
dir.create(dirname(fig_path), recursive = TRUE, showWarnings = FALSE)
ggsave(fig_path, p, width = 10, height = 5, dpi = 150)
cat(sprintf("Figure saved: %s\n", fig_path))

# Save summary
saveRDS(summary_df, file.path(cfg$paths$calibration_dir, "calibration_summary.rds"))
write.csv(summary_df, file.path(cfg$paths$calibration_dir, "calibration_summary.csv"),
          row.names = FALSE)

cat("\nPhase 1a COMPLETE.\n")
