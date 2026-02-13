#!/usr/bin/env Rscript
# =============================================================================
# test_pilot.R  —  Pilot Test with Parallel Computing + IRW Items
# =============================================================================
# Runs 4 strategically selected conditions with full MCMC chains (10K iter)
# to verify: (1) all code works end-to-end, (2) results match theory,
# (3) timing estimates for full simulation, (4) parallel computing works.
#
# Parallelization strategy: condition-SEQUENTIAL, replication-PARALLEL
#   - Outer loop: iterate over conditions one at a time (sequential)
#   - Inner loop: distribute K replications across workers (parallel)
#   - Each worker independently compiles NIMBLE model (no shared state)
#   - Uses future::multisession to spawn separate R processes
#
# Item source: IRW (Item Response Warehouse) — set in sim_config.yaml
#
# Selected conditions (covering key variation):
#   C1: Rasch, normal, N=100, rho=0.70  (baseline)
#   C2: Rasch, bimodal, N=200, rho=0.80 (non-normal, larger N)
#   C3: 2PL,   normal, N=100, rho=0.50  (low reliability)
#   C4: 2PL,   skew_pos, N=50, rho=0.90 (high reliability, small N, non-normal)
#
# Total fits: 4 conditions x 3 reps x 3 priors = 36 MCMC fits
# =============================================================================

cat("\n====================================================\n")
cat("  PILOT TEST: 4 conditions x 3 reps x 3 priors\n")
cat("  Full MCMC (10K iter, 2K burnin)\n")
cat("  WITH PARALLEL COMPUTING + IRW ITEMS\n")
cat("====================================================\n\n")

library(nimble)
library(DPMirt)
library(DPprior)
library(IRTsimrel)
library(tidyverse)
library(yaml)
library(digest)
library(future)
library(future.apply)

source("R/00_config.R")
source("R/01_calibrate.R")
source("R/02_generate.R")
source("R/03_fit.R")
source("R/04_evaluate.R")

cfg <- load_config()

# Use 3 reps for pilot (not 100)
K_pilot <- 3L
cfg$posterior_storage$save_posteriors <- FALSE

# Verify IRW item source
cat(sprintf("Item source: %s\n", cfg$calibration$item_source))
stopifnot(cfg$calibration$item_source == "irw")

# Clean old checkpoints
if (dir.exists(cfg$paths$checkpoint_dir)) {
  old_ckpts <- list.files(cfg$paths$checkpoint_dir, full.names = TRUE)
  if (length(old_ckpts) > 0) file.remove(old_ckpts)
}

# =============================================================================
# Setup
# =============================================================================
design <- create_design_matrix(cfg)

# Select 4 pilot conditions
pilot_conditions <- design %>%
  filter(
    (irt_model == "rasch" & latent_shape == "normal"   & N == 100 & target_rho == 0.70) |
    (irt_model == "rasch" & latent_shape == "bimodal"  & N == 200 & target_rho == 0.80) |
    (irt_model == "2pl"   & latent_shape == "normal"   & N == 100 & target_rho == 0.50) |
    (irt_model == "2pl"   & latent_shape == "skew_pos" & N == 50  & target_rho == 0.90)
  )

cat(sprintf("\nPilot conditions: %d\n", nrow(pilot_conditions)))
for (i in 1:nrow(pilot_conditions)) {
  r <- pilot_conditions[i, ]
  cat(sprintf("  C%d: %s, %s, N=%d, rho=%.2f (id=%d)\n",
      i, r$irt_model, r$latent_shape, r$N, r$target_rho, r$condition_id))
}
cat(sprintf("\nTotal fits: %d conditions x %d reps x 3 priors = %d\n\n",
    nrow(pilot_conditions), K_pilot, nrow(pilot_conditions) * K_pilot * 3))

# =============================================================================
# Phase 1a: Calibrate all cells (with IRW items)
# =============================================================================
cat("=== Phase 1a: EQC Calibration (IRW item source) ===\n")
calib_results <- run_all_calibrations(cfg)

# Quick verification of IRW items
cat("\n--- IRW Item Parameter Checks ---\n")
sample_cell <- calib_results[[1]]
cat(sprintf("  Sample cell: %s\n", names(calib_results)[1]))
cat(sprintf("  beta range:   [%.3f, %.3f], mean=%.3f, sd=%.3f\n",
    min(sample_cell$beta_vec), max(sample_cell$beta_vec),
    mean(sample_cell$beta_vec), sd(sample_cell$beta_vec)))
cat(sprintf("  lambda range: [%.3f, %.3f], mean=%.3f, sd=%.3f\n",
    min(sample_cell$lambda_scaled), max(sample_cell$lambda_scaled),
    mean(sample_cell$lambda_scaled), sd(sample_cell$lambda_scaled)))

# Build alpha prior cache
alpha_cache <- build_alpha_cache(cfg)
cat(sprintf("\nAlpha cache: %d entries\n", length(alpha_cache)))

# =============================================================================
# Setup parallel workers
# =============================================================================
n_workers <- cfg$parallel$n_workers
cat(sprintf("\n=== Setting up parallel computing: %d workers (%s) ===\n",
    n_workers, cfg$parallel$strategy))

plan(multisession, workers = n_workers)
cat(sprintf("Workers available: %d\n", nbrOfWorkers()))
cat("Parallel backend ready.\n\n")

# =============================================================================
# Run pilot: SEQUENTIAL over conditions, PARALLEL over replications
# =============================================================================
all_pilot_results <- list()
timing_log <- list()

t_global_start <- Sys.time()

# Get project root and source file paths for workers
project_root <- cfg$paths$project_root

for (ci in 1:nrow(pilot_conditions)) {
  row <- pilot_conditions[ci, ]
  eqc <- lookup_calibration(row, calib_results)

  cat(sprintf("=== Condition %d/%d: %s, %s, N=%d, rho=%.2f ===\n",
      ci, nrow(pilot_conditions),
      row$irt_model, row$latent_shape, row$N, row$target_rho))

  # --- Parallel execution of K replications ---
  t_cond_start <- Sys.time()

  rep_results <- future_lapply(1:K_pilot, function(k) {
    # Each worker needs to load packages and source files independently
    # because multisession spawns new R processes
    suppressMessages({
      library(nimble)
      library(DPMirt)
      library(DPprior)
      library(IRTsimrel)
      library(tidyverse)
      library(yaml)
      library(digest)
    })

    setwd(project_root)
    source("R/00_config.R")
    source("R/01_calibrate.R")
    source("R/02_generate.R")
    source("R/03_fit.R")
    source("R/04_evaluate.R")

    # Load config inside worker
    w_cfg <- load_config()
    w_cfg$posterior_storage$save_posteriors <- FALSE

    t0 <- Sys.time()

    result <- run_one_replication(
      design_row  = row,
      eqc_result  = eqc,
      rep_id      = k,
      alpha_cache = alpha_cache,
      config      = w_cfg
    )

    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    list(result = result, time_secs = elapsed, rep_id = k)
  }, future.seed = TRUE)

  t_cond_elapsed <- as.numeric(difftime(Sys.time(), t_cond_start, units = "secs"))

  # Collect results from workers
  for (res in rep_results) {
    cat(sprintf("  Rep %d: %.1f sec (all 3 priors)\n", res$rep_id, res$time_secs))

    all_pilot_results[[length(all_pilot_results) + 1]] <- res$result
    timing_log[[length(timing_log) + 1]] <- data.frame(
      condition_id = row$condition_id,
      irt_model = row$irt_model,
      N = row$N,
      rep_id = res$rep_id,
      time_secs = res$time_secs
    )
  }

  cat(sprintf("  Wall-clock for condition: %.1f sec (%d reps in parallel)\n\n",
      t_cond_elapsed, K_pilot))
}

t_global_end <- Sys.time()
total_time <- as.numeric(difftime(t_global_end, t_global_start, units = "mins"))

# Shut down parallel workers
plan(sequential)

# =============================================================================
# Compile results
# =============================================================================
cat("====================================================\n")
cat("  COMPILING RESULTS\n")
cat("====================================================\n\n")

results_df <- dplyr::bind_rows(all_pilot_results)
timing_df <- dplyr::bind_rows(timing_log)

# Focus on theta losses
theta_results <- results_df %>%
  filter(parameter == "theta") %>%
  mutate(
    rmse = sqrt(msel),
    cond_label = sprintf("%s_%s_N%d_r%02d",
                         irt_model, latent_shape, N, as.integer(target_rho * 100))
  )

# Summary by condition x prior x method
summary_df <- theta_results %>%
  group_by(irt_model, latent_shape, N, target_rho, prior, method) %>%
  summarise(
    mean_rmse  = mean(rmse, na.rm = TRUE),
    mean_mselr = mean(mselr, na.rm = TRUE),
    mean_ks    = mean(ks, na.rm = TRUE),
    n_reps     = n(),
    .groups = "drop"
  )

cat("=== THETA LOSS SUMMARY (averaged over reps) ===\n\n")
for (ci in 1:nrow(pilot_conditions)) {
  r <- pilot_conditions[ci, ]
  cat(sprintf("Condition: %s, %s, N=%d, rho=%.2f\n",
      r$irt_model, r$latent_shape, r$N, r$target_rho))
  cat(sprintf("%-12s %-4s  %8s  %10s  %8s\n", "Prior", "Meth", "RMSE", "MSELR", "KS"))
  cat(paste(rep("-", 55), collapse = ""), "\n")

  cond_summary <- summary_df %>%
    filter(irt_model == r$irt_model, latent_shape == r$latent_shape,
           N == r$N, target_rho == r$target_rho)

  for (j in 1:nrow(cond_summary)) {
    s <- cond_summary[j, ]
    cat(sprintf("%-12s %-4s  %8.4f  %10.6f  %8.4f\n",
        s$prior, toupper(s$method), s$mean_rmse, s$mean_mselr, s$mean_ks))
  }
  cat("\n")
}

# =============================================================================
# Timing estimates for full simulation
# =============================================================================
cat("====================================================\n")
cat("  TIMING ESTIMATES\n")
cat("====================================================\n\n")

mean_time_per_rep <- mean(timing_df$time_secs)
cat(sprintf("Mean time per replication (3 priors): %.1f sec\n", mean_time_per_rep))
cat(sprintf("Actual pilot wall-clock: %.1f min\n", total_time))

# Compute parallel speedup
sequential_estimate <- sum(timing_df$time_secs) / 60
speedup <- sequential_estimate / total_time
cat(sprintf("Sequential equivalent:   %.1f min\n", sequential_estimate))
cat(sprintf("Parallel speedup:        %.2fx (%d workers)\n", speedup, n_workers))

time_by_N <- timing_df %>%
  group_by(N) %>%
  summarise(mean_secs = mean(time_secs), .groups = "drop")
cat("\nTime by N (per replication):\n")
for (i in 1:nrow(time_by_N)) {
  cat(sprintf("  N=%3d: %.1f sec/rep\n", time_by_N$N[i], time_by_N$mean_secs[i]))
}

# Estimate full simulation
total_conditions <- 120
total_reps <- cfg$design$replications
estimated_hours_seq <- (mean_time_per_rep * total_conditions * total_reps) / 3600
cat(sprintf("\nFull simulation estimates:\n"))
cat(sprintf("  Sequential: %d conditions x %d reps x %.1f sec/rep = %.0f hours\n",
    total_conditions, total_reps, mean_time_per_rep, estimated_hours_seq))
cat(sprintf("  With %d parallel workers: ~%.0f hours\n",
    n_workers, estimated_hours_seq / speedup))

cat(sprintf("\nTotal pilot wall-clock time: %.1f minutes\n", total_time))

# =============================================================================
# Save results
# =============================================================================
saveRDS(results_df, file.path(cfg$paths$results_dir, "pilot_results_raw.rds"))
saveRDS(summary_df, file.path(cfg$paths$results_dir, "pilot_results_summary.rds"))
write.csv(summary_df, file.path(cfg$paths$results_dir, "pilot_results_summary.csv"),
          row.names = FALSE)

# =============================================================================
# Sanity checks
# =============================================================================
cat("\n====================================================\n")
cat("  SANITY CHECKS\n")
cat("====================================================\n\n")

# Check 1: All fits succeeded
n_failed <- sum(theta_results$status != "success", na.rm = TRUE)
cat(sprintf("[%s] All fits succeeded (%d/%d)\n",
    ifelse(n_failed == 0, "PASS", "FAIL"),
    sum(theta_results$status == "success"), nrow(theta_results)))

# Check 2: Higher reliability -> lower RMSE (for PM, Gaussian, normal)
gauss_pm <- summary_df %>%
  filter(prior == "gaussian", method == "pm",
         irt_model == "rasch", latent_shape == "normal")
if (nrow(gauss_pm) > 0) {
  rmse_decreasing <- all(diff(gauss_pm$mean_rmse[order(gauss_pm$target_rho)]) < 0)
  cat(sprintf("[%s] RMSE decreases with reliability (Gauss PM, Rasch, Normal)\n",
      ifelse(rmse_decreasing, "PASS", "NOTE")))
}

# Check 3: GR has best KS across conditions (on average)
ks_by_method <- theta_results %>%
  group_by(method) %>%
  summarise(mean_ks = mean(ks, na.rm = TRUE), .groups = "drop")
best_ks_method <- ks_by_method$method[which.min(ks_by_method$mean_ks)]
cat(sprintf("[%s] Best KS method overall: %s (theory predicts GR)\n",
    ifelse(best_ks_method == "gr", "PASS", "NOTE"), toupper(best_ks_method)))

# Check 4: PM has lowest MSEL (it's Bayes optimal for MSE)
msel_by_method <- theta_results %>%
  group_by(method) %>%
  summarise(mean_msel = mean(msel, na.rm = TRUE), .groups = "drop")
best_msel_method <- msel_by_method$method[which.min(msel_by_method$mean_msel)]
cat(sprintf("[%s] Best MSEL method overall: %s (theory predicts PM)\n",
    ifelse(best_msel_method == "pm", "PASS", "NOTE"), toupper(best_msel_method)))

# Check 5: Parallel speedup is reasonable (> 1.5x with multi-workers)
cat(sprintf("[%s] Parallel speedup: %.2fx (expected > 1.5x with %d workers)\n",
    ifelse(speedup > 1.5, "PASS", "NOTE"), speedup, n_workers))

cat("\n====================================================\n")
cat("  PILOT TEST COMPLETE\n")
cat(sprintf("  Item source: %s\n", cfg$calibration$item_source))
cat(sprintf("  Parallel workers: %d\n", n_workers))
cat(sprintf("  Parallel speedup: %.2fx\n", speedup))
cat(sprintf("  Wall-clock time: %.1f minutes\n", total_time))
cat("====================================================\n")
