#!/usr/bin/env Rscript
# =============================================================================
# run_simulation.R  —  Master Orchestration Script (Parallel)
# =============================================================================
# Runs the full simulation study with configurable scope:
#   --pilot      : 4 conditions x 3 reps (for testing)
#   --subset N   : First N conditions x full reps
#   (default)    : All 120 conditions x K reps
#
# Parallelization: condition-SEQUENTIAL, replication-PARALLEL
#   Uses future::multisession to distribute replications across workers.
#   Each worker independently loads packages and compiles NIMBLE models.
#
# Usage:
#   Rscript scripts/run_simulation.R --pilot
#   Rscript scripts/run_simulation.R --subset 10
#   Rscript scripts/run_simulation.R
#
# The script is checkpoint-aware: it skips replications that already have
# checkpoint files on disk.
# =============================================================================

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
source("R/05_postprocess.R")

# =============================================================================
# Parse command line arguments
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)

pilot_mode <- "--pilot" %in% args
subset_n <- NA
if ("--subset" %in% args) {
  idx <- which(args == "--subset")
  if (idx < length(args)) subset_n <- as.integer(args[idx + 1])
}

sequential_mode <- "--sequential" %in% args

# =============================================================================
# Load config and setup
# =============================================================================
cfg <- load_config()
ensure_directories(cfg)

if (pilot_mode) {
  K <- 3L
  cfg$posterior_storage$save_posteriors <- FALSE
  cat("=== PILOT MODE: 4 conditions x 3 reps ===\n\n")
} else {
  K <- cfg$design$replications
}

design <- create_design_matrix(cfg)

# Select conditions
if (pilot_mode) {
  conditions <- design %>%
    filter(
      (irt_model == "rasch" & latent_shape == "normal"   & N == 100 & target_rho == 0.70) |
      (irt_model == "rasch" & latent_shape == "bimodal"  & N == 200 & target_rho == 0.80) |
      (irt_model == "2pl"   & latent_shape == "normal"   & N == 100 & target_rho == 0.50) |
      (irt_model == "2pl"   & latent_shape == "skew_pos" & N == 50  & target_rho == 0.90)
    )
} else if (!is.na(subset_n)) {
  conditions <- design %>% slice(1:min(subset_n, nrow(design)))
} else {
  conditions <- design
}

n_cond <- nrow(conditions)
total_fits <- n_cond * K * 3L

cat(sprintf("Conditions: %d\n", n_cond))
cat(sprintf("Replications per condition: %d\n", K))
cat(sprintf("Total MCMC fits: %d\n", total_fits))
cat(sprintf("Item source: %s\n\n", cfg$calibration$item_source))

# =============================================================================
# Phase 1a: Calibration
# =============================================================================
cat("=== Phase 1a: EQC Calibration ===\n")
calib_results <- run_all_calibrations(cfg)

# =============================================================================
# Pre-compute alpha prior cache
# =============================================================================
cat("\n=== Pre-computing alpha prior cache ===\n")
alpha_cache <- build_alpha_cache(cfg)
cat(sprintf("Cached %d DP prior specifications\n\n", length(alpha_cache)))

# =============================================================================
# Setup parallel workers
# =============================================================================
project_root <- cfg$paths$project_root

if (!sequential_mode) {
  n_workers <- cfg$parallel$n_workers
  cat(sprintf("=== Parallel computing: %d workers (%s) ===\n",
      n_workers, cfg$parallel$strategy))
  plan(multisession, workers = n_workers)
  cat(sprintf("Workers active: %d\n\n", nbrOfWorkers()))
} else {
  cat("=== Sequential mode (--sequential flag) ===\n\n")
}

# =============================================================================
# Main loop: condition x replication (parallel within condition)
# =============================================================================
cat("=== Main Simulation Loop ===\n\n")

t_global <- Sys.time()
fit_count <- 0L
skip_count <- 0L

for (ci in seq_len(n_cond)) {
  row <- conditions[ci, ]
  eqc <- lookup_calibration(row, calib_results)

  # Check which reps are already done
  done_reps <- get_completed_reps(row$condition_id, cfg)

  cat(sprintf("[%d/%d] %s %s N=%d rho=%.2f (id=%d)",
      ci, n_cond, row$irt_model, row$latent_shape, row$N,
      row$target_rho, row$condition_id))

  if (length(done_reps) >= K) {
    cat(sprintf(" -- SKIP (all %d reps done)\n", K))
    skip_count <- skip_count + K * 3L
    next
  }

  remaining <- setdiff(1:K, done_reps)
  cat(sprintf(" -- %d reps remaining\n", length(remaining)))

  if (!sequential_mode && length(remaining) > 1) {
    # --- PARALLEL execution ---
    t_cond <- Sys.time()

    rep_results <- future_lapply(remaining, function(k) {
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

      w_cfg <- load_config()
      if (pilot_mode) w_cfg$posterior_storage$save_posteriors <- FALSE

      t0 <- Sys.time()

      result <- tryCatch(
        run_one_replication(
          design_row  = row,
          eqc_result  = eqc,
          rep_id      = k,
          alpha_cache = alpha_cache,
          config      = w_cfg
        ),
        error = function(e) {
          warning(sprintf("  ERROR in cond %d rep %d: %s",
                          row$condition_id, k, conditionMessage(e)))
          NULL
        }
      )

      elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      list(result = result, elapsed = elapsed, rep_id = k)
    }, future.seed = TRUE)

    t_cond_elapsed <- as.numeric(difftime(Sys.time(), t_cond, units = "secs"))

    for (res in rep_results) {
      n_failed <- if (!is.null(res$result)) sum(res$result$status != "success") else 3L
      cat(sprintf("    rep %d: %.1f sec%s\n", res$rep_id, res$elapsed,
          if (n_failed > 0) sprintf(" (%d FAILED)", n_failed) else ""))
      fit_count <- fit_count + 3L
    }
    cat(sprintf("    Wall-clock: %.1f sec (%d reps parallel)\n",
        t_cond_elapsed, length(remaining)))

  } else {
    # --- SEQUENTIAL execution ---
    for (k in remaining) {
      t0 <- Sys.time()

      result <- tryCatch(
        run_one_replication(
          design_row  = row,
          eqc_result  = eqc,
          rep_id      = k,
          alpha_cache = alpha_cache,
          config      = cfg
        ),
        error = function(e) {
          warning(sprintf("  ERROR in cond %d rep %d: %s",
                          row$condition_id, k, conditionMessage(e)))
          NULL
        }
      )

      elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      fit_count <- fit_count + 3L

      n_failed <- if (!is.null(result)) sum(result$status != "success") else 3L
      cat(sprintf("    rep %d: %.1f sec%s\n", k, elapsed,
          if (n_failed > 0) sprintf(" (%d FAILED)", n_failed) else ""))

      # Periodic GC
      if (fit_count %% 30 == 0) gc(verbose = FALSE)
    }
  }
}

t_total <- as.numeric(difftime(Sys.time(), t_global, units = "mins"))

# Shut down parallel workers
if (!sequential_mode) plan(sequential)

cat(sprintf("\n=== Simulation Complete ===\n"))
cat(sprintf("Fits completed: %d\n", fit_count))
cat(sprintf("Fits skipped:   %d\n", skip_count))
cat(sprintf("Total time:     %.1f minutes\n\n", t_total))

# =============================================================================
# Compile results
# =============================================================================
cat("=== Compiling Results ===\n")
all_results <- compile_all_results(cfg)
cat(sprintf("Total result rows: %d\n", nrow(all_results)))

# Save compiled results
saveRDS(all_results, file.path(cfg$paths$results_dir, "simulation_results.rds"))
write.csv(all_results, file.path(cfg$paths$results_dir, "simulation_results.csv"),
          row.names = FALSE)
cat("Results saved to output/results/\n")

cat("\nDone.\n")
