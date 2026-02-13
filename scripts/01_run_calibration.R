# =============================================================================
# 01_run_calibration.R
# =============================================================================
# Phase 1: Run EQC calibration for all conditions in the design matrix.
#
# This script runs sequentially (calibration is fast and memory-light).
# Each calibration result is cached to disk; re-running this script will
# skip already-calibrated conditions.
#
# Prerequisites: Run 00_setup.R first (or source all R modules).
#
# Author: JoonHo Lee
# Date: February 2026
# =============================================================================


# --- Source setup if not already loaded ---
if (!exists("config")) {
  source(file.path(here::here(), "scripts", "00_setup.R"))
}


# =============================================================================
# Run Calibration
# =============================================================================

log_msg("Starting calibration phase...", config = config)

tic <- proc.time()

calib_results <- run_all_calibrations(design, config)

elapsed <- (proc.time() - tic)[["elapsed"]]
log_msg("Calibration phase complete in ", format_time(elapsed), config = config)


# =============================================================================
# Verify Calibration Quality
# =============================================================================

log_msg("Verifying calibration quality...", config = config)

# Build a summary table of calibration results
calib_summary <- tibble::tibble(
  condition_id = design$condition_id,
  irt_model    = design$irt_model,
  latent_shape = design$latent_shape,
  N            = design$N,
  target_rho   = design$target_rho,
  c_star       = vapply(calib_results, function(cr) cr$c_star, numeric(1)),
  achieved_rho = vapply(calib_results, function(cr) cr$achieved_rho, numeric(1)),
  rho_error    = vapply(calib_results, function(cr) {
    abs(cr$achieved_rho - cr$target_rho)
  }, numeric(1))
)

# Report calibration accuracy
cat("\n--- Calibration Summary ---\n")
cat(sprintf("  Total conditions: %d\n", nrow(calib_summary)))
cat(sprintf("  Max |rho_error|:  %.6f\n", max(calib_summary$rho_error)))
cat(sprintf("  Mean |rho_error|: %.6f\n", mean(calib_summary$rho_error)))

# Flag any poor calibrations
poor_calib <- calib_summary[calib_summary$rho_error > 0.01, ]
if (nrow(poor_calib) > 0) {
  warning(nrow(poor_calib), " conditions with |rho_error| > 0.01:")
  print(poor_calib)
} else {
  cat("  All calibrations within tolerance.\n")
}

# Save summary
calib_summary_path <- file.path(
  config$paths$calibration_dir, "calibration_summary.rds"
)
saveRDS(calib_summary, calib_summary_path)
log_msg("Calibration summary saved to: ", calib_summary_path, config = config)

cat("\nCalibration phase complete. Proceed to 02_run_simulation.R\n")
