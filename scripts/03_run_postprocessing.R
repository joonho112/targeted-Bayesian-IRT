# =============================================================================
# 03_run_postprocessing.R
# =============================================================================
# Phase 3: Post-processing of simulation results.
#
# Steps:
#   1. Load compiled results
#   2. Compute summary statistics (mean/median/SD of losses per condition)
#   3. Compute scaled EDF (pooled across replications)
#   4. Fit meta-model regressions
#   5. Extract marginal effects
#   6. Save all post-processed outputs
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
# 3. Compute Scaled EDF
# =============================================================================

log_msg("Computing scaled EDFs (this may take a while)...", config = config)

tic_edf <- proc.time()

edf_results <- compute_scaled_edf(
  results_df    = results_all,
  posteriors_dir = config$paths$posteriors_dir,
  config        = config
)

edf_elapsed <- (proc.time() - tic_edf)[["elapsed"]]
log_msg("Scaled EDF computation complete in ", format_time(edf_elapsed),
        config = config)

# Save
edf_path <- file.path(config$paths$results_dir, "scaled_edf.rds")
saveRDS(edf_results, edf_path)
log_msg("Scaled EDF saved to: ", edf_path, config = config)


# =============================================================================
# 4. Run Meta-Model Regressions
# =============================================================================

log_msg("Fitting meta-model regressions...", config = config)

tic_reg <- proc.time()

meta_reg_results <- run_meta_regression(results_all, config)

reg_elapsed <- (proc.time() - tic_reg)[["elapsed"]]
log_msg("Meta-regression complete in ", format_time(reg_elapsed),
        config = config)

# Save
reg_path <- file.path(config$paths$results_dir, "meta_regression_results.rds")
saveRDS(meta_reg_results, reg_path)
log_msg("Meta-regression results saved to: ", reg_path, config = config)


# =============================================================================
# 5. Extract Marginal Effects
# =============================================================================

log_msg("Extracting marginal effects...", config = config)

marginal_effects <- extract_marginal_effects(meta_reg_results)

# Save
me_path <- file.path(config$paths$results_dir, "marginal_effects.rds")
saveRDS(marginal_effects, me_path)
log_msg("Marginal effects saved to: ", me_path, config = config)


# =============================================================================
# 6. Print Summary of Meta-Regression Results
# =============================================================================

cat("\n--- Meta-Regression Summary ---\n\n")

for (i in seq_len(nrow(meta_reg_results))) {
  row <- meta_reg_results[i, ]
  cat(sprintf(
    "  [%s | %s | %s] R^2 = %.3f (n = %d)\n",
    row$metric, row$latent_shape, row$irt_model,
    row$r_squared, row$n_obs
  ))
}

cat("\n--- Marginal Effects (top 10 by |estimate|) ---\n\n")
me_sorted <- marginal_effects[order(-abs(marginal_effects$estimate)), ]
print(head(me_sorted, 10), n = 10)


# =============================================================================
# 7. Summary
# =============================================================================

cat("\n=== Post-Processing Complete ===\n\n")
cat("Outputs saved:\n")
cat(sprintf("  Summary (mean):     %s\n",
            file.path(config$paths$results_dir, "results_summary_mean.rds")))
cat(sprintf("  Summary (trimmed):  %s\n",
            file.path(config$paths$results_dir, "results_summary_trimmed.rds")))
cat(sprintf("  Scaled EDF:         %s\n", edf_path))
cat(sprintf("  Meta-regression:    %s\n", reg_path))
cat(sprintf("  Marginal effects:   %s\n", me_path))
cat("\nProceed to 04_generate_figures.R\n")
