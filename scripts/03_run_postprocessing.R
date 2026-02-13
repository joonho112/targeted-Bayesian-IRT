# =============================================================================
# 03_run_postprocessing.R
# =============================================================================
# Phase 3: Post-processing of simulation results.
#
# Steps:
#   1. Load compiled results
#   2. Compute summary statistics (mean/median/SD of losses per condition)
#   3. Fit meta-model regressions (per metric)
#   4. Save all post-processed outputs
#
# Prerequisites:
#   - 02_run_simulation.R has been run (compiled results on disk).
#
# Author: JoonHo Lee
# Date: February 2026
# =============================================================================


# --- Source setup if not already loaded ---
if (!exists("config")) {
  source(file.path(here::here(), "scripts", "00_setup.R"))
}


# =============================================================================
# 1. Load Compiled Results
# =============================================================================

results_path <- file.path(
  config$paths$results_dir, "simulation_results_compiled.rds"
)

if (!file.exists(results_path)) {
  log_msg("Compiled results not found. Compiling from checkpoints...",
          config = config)
  results_all <- compile_all_results(config)
} else {
  log_msg("Loading compiled results from: ", results_path, config = config)
  results_all <- readRDS(results_path)
}

cat(sprintf("Loaded %d result rows.\n", nrow(results_all)))


# =============================================================================
# 2. Compute Summary Statistics
# =============================================================================

log_msg("Computing summary statistics...", config = config)

# Untrimmed mean
results_summary <- summarize_results(results_all, trim = 0)

# 10% trimmed mean (robust to outliers)
results_summary_trimmed <- summarize_results(results_all, trim = 0.10)

# Save summaries
saveRDS(results_summary,
        file.path(config$paths$results_dir, "results_summary_mean.rds"))
saveRDS(results_summary_trimmed,
        file.path(config$paths$results_dir, "results_summary_trimmed.rds"))

log_msg("Summary statistics saved.", config = config)


# =============================================================================
# 3. Run Meta-Model Regressions
# =============================================================================

log_msg("Fitting meta-model regressions...", config = config)

tic_reg <- proc.time()

meta_results <- list()
for (metric in config$metrics) {
  log_msg("  Fitting meta-model for: ", metric, config = config)
  meta_results[[metric]] <- run_meta_regression(
    results_df = results_all,
    metric     = metric,
    config     = config
  )
}

reg_elapsed <- (proc.time() - tic_reg)[["elapsed"]]
log_msg("Meta-regression complete in ", format_time(reg_elapsed),
        config = config)

# Save
reg_path <- file.path(config$paths$results_dir, "meta_regression_results.rds")
saveRDS(meta_results, reg_path)
log_msg("Meta-regression results saved to: ", reg_path, config = config)


# =============================================================================
# 4. Print Summary of Meta-Regression Results
# =============================================================================

cat("\n--- Meta-Regression Summary ---\n\n")

for (metric in names(meta_results)) {
  res <- meta_results[[metric]]
  if (!is.null(res)) {
    cat(sprintf(
      "  [%s] R^2 = %.3f, Adj-R^2 = %.3f (n = %d)\n",
      metric, res$r_squared, res$adj_r_sq, res$n_obs
    ))
  } else {
    cat(sprintf("  [%s] Skipped (insufficient data)\n", metric))
  }
}


# =============================================================================
# 5. Summary
# =============================================================================

cat("\n=== Post-Processing Complete ===\n\n")
cat("Outputs saved:\n")
cat(sprintf("  Summary (mean):     %s\n",
            file.path(config$paths$results_dir, "results_summary_mean.rds")))
cat(sprintf("  Summary (trimmed):  %s\n",
            file.path(config$paths$results_dir, "results_summary_trimmed.rds")))
cat(sprintf("  Meta-regression:    %s\n", reg_path))
cat("\nProceed to 04_generate_figures.R\n")
