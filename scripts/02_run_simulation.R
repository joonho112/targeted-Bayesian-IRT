# =============================================================================
# 02_run_simulation.R
# =============================================================================
# Phase 2: Main simulation loop.
#
# Strategy:
#   - Process CONDITIONS sequentially (memory management: one condition's
#     calibration result + datasets in memory at a time).
#   - Within each condition, REPLICATIONS are parallelized using
#     future + future.apply (the unit of parallelization is run_one_replication).
#   - Checkpointing: After each replication, results are saved to an
#     incremental checkpoint file. On resume, completed replications
#     are detected and skipped.
#
# Prerequisites:
#   - 00_setup.R has been run (packages loaded, modules sourced).
#   - 01_run_calibration.R has been run (calibration results on disk).
#
# Author: JoonHo Lee
# Date: February 2026
# =============================================================================


# --- Source setup if not already loaded ---
if (!exists("config")) {
  source(file.path(here::here(), "scripts", "00_setup.R"))
}


# =============================================================================
# 1. Load Calibration Results & Build Design
# =============================================================================

log_msg("Loading calibration results...", config = config)

calib_results <- run_all_calibrations(config)   # uses cache; returns named list

design <- create_design_matrix(config)

log_msg("Loaded ", length(calib_results), " calibration cells for ",
        nrow(design), " design conditions.", config = config)


# =============================================================================
# 2. Pre-compute Alpha Cache (DP prior hyperparameters)
# =============================================================================

log_msg("Building alpha prior cache...", config = config)
alpha_cache <- build_alpha_cache(config)
log_msg("Alpha cache ready: ", length(alpha_cache), " entries.", config = config)


# =============================================================================
# 3. Configure Parallel Backend
# =============================================================================

n_workers <- config$parallel$n_workers
log_msg("Setting up parallel backend: ", config$parallel$strategy,
        " with ", n_workers, " workers.", config = config)

future::plan(future::multisession, workers = n_workers)

project_root <- here::here()


# =============================================================================
# 3.5  Quick Parallel Smoke Test (first condition, n_workers reps)
# =============================================================================
# Before committing to the full loop, verify that parallel workers can
# independently load packages, compile nimble models, and return results.

cat("\n")
log_msg("=== Parallel Smoke Test ===", config = config)

test_row   <- design[1, ]
test_eqc   <- lookup_calibration(test_row, calib_results)
test_K     <- n_workers                          # one rep per worker
test_reps  <- seq_len(test_K)

log_msg("  Condition: ", test_row$condition_id,
        " (", test_row$irt_model, "/", test_row$latent_shape,
        "/N=", test_row$N, "/rho=", test_row$target_rho, ")",
        config = config)
log_msg("  Launching ", test_K, " replications on ", n_workers,
        " workers ...", config = config)

smoke_tic <- proc.time()

smoke_results <- future.apply::future_lapply(
  test_reps,
  function(rep_id) {
    suppressMessages({
      library(nimble); library(DPMirt); library(DPprior)
      library(IRTsimrel); library(tidyverse); library(yaml); library(digest)
    })
    setwd(project_root)
    source("R/utils.R");    source("R/00_config.R")
    source("R/01_calibrate.R"); source("R/02_generate.R")
    source("R/03_fit.R");   source("R/04_evaluate.R")

    w_cfg  <- load_config()
    t0     <- proc.time()

    result <- run_one_replication(
      design_row  = test_row,
      eqc_result  = test_eqc,
      rep_id      = rep_id,
      alpha_cache = alpha_cache,
      config      = w_cfg
    )

    elapsed <- (proc.time() - t0)[["elapsed"]]
    list(result = result, time_secs = elapsed, rep_id = rep_id)
  },
  future.seed = TRUE
)

smoke_elapsed <- (proc.time() - smoke_tic)[["elapsed"]]
smoke_times   <- vapply(smoke_results, function(x) x$time_secs, numeric(1))
smoke_ok      <- vapply(smoke_results, function(x) {
  all(x$result$status == "success")
}, logical(1))

log_msg("  Wall-clock:  ", format_time(smoke_elapsed), config = config)
log_msg("  Per-worker:  ", sprintf("%.1f–%.1fs", min(smoke_times), max(smoke_times)),
        config = config)
log_msg("  Success:     ", sum(smoke_ok), "/", test_K, config = config)

if (!all(smoke_ok)) {
  stop("Parallel smoke test FAILED. Fix errors before running full simulation.")
}
log_msg("=== Smoke Test PASSED ===", config = config)
cat("\n")

# Clean up smoke-test checkpoints so they don't pollute the main run
smoke_ckpt_dir <- file.path(config$paths$results_dir, "checkpoints")
smoke_files <- list.files(
  smoke_ckpt_dir,
  pattern = sprintf("^checkpoint_%04d_rep", test_row$condition_id),
  full.names = TRUE
)
if (length(smoke_files) > 0) file.remove(smoke_files)


# =============================================================================
# 4. Main Simulation Loop (Conditions x Replications)
# =============================================================================

K <- config$design$replications
n_conditions <- nrow(design)

log_msg(
  "Starting simulation: ", n_conditions, " conditions x ",
  K, " replications x 3 priors = ",
  n_conditions * K * 3, " total model fits.",
  config = config
)

overall_tic <- proc.time()

for (ci in seq_len(n_conditions)) {

  design_row <- design[ci, ]
  cond_id    <- design_row$condition_id

  # Look up the correct EQC calibration result for this condition
  eqc_result <- lookup_calibration(design_row, calib_results)

  # --- Check which replications are already done ---
  completed_reps <- get_completed_reps(cond_id, config)
  remaining_reps <- setdiff(seq_len(K), completed_reps)

  if (length(remaining_reps) == 0) {
    log_msg("Condition ", cond_id, ": all ", K,
            " replications already complete. Skipping.", config = config)
    next
  }

  log_msg(
    "Condition ", cond_id, " (", ci, "/", n_conditions, "): ",
    length(remaining_reps), " of ", K, " replications remaining. ",
    "(model=", design_row$irt_model,
    ", shape=", design_row$latent_shape,
    ", N=", design_row$N,
    ", rho=", design_row$target_rho, ")",
    config = config
  )

  cond_tic <- proc.time()

  # --- Parallelize replications within this condition ---
  rep_results <- future.apply::future_lapply(
    remaining_reps,
    function(rep_id) {

      # Each worker must independently load packages and source modules
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
      source("R/utils.R")
      source("R/00_config.R")
      source("R/01_calibrate.R")
      source("R/02_generate.R")
      source("R/03_fit.R")
      source("R/04_evaluate.R")

      w_cfg <- load_config()

      # Run one complete replication
      result <- run_one_replication(
        design_row  = design_row,
        eqc_result  = eqc_result,
        rep_id      = rep_id,
        alpha_cache = alpha_cache,
        config      = w_cfg
      )

      result
    },
    future.seed = TRUE
  )

  cond_elapsed <- (proc.time() - cond_tic)[["elapsed"]]
  log_msg(
    "Condition ", cond_id, " complete in ", format_time(cond_elapsed),
    " (", length(remaining_reps), " reps).",
    config = config
  )
}

overall_elapsed <- (proc.time() - overall_tic)[["elapsed"]]
log_msg("Simulation phase complete in ", format_time(overall_elapsed),
        config = config)


# =============================================================================
# 5. Compile All Results
# =============================================================================

log_msg("Compiling all checkpoint files...", config = config)

results_all <- compile_all_results(config)

cat(sprintf("\n--- Simulation Complete ---\n"))
cat(sprintf("  Total rows:     %d\n", nrow(results_all)))
cat(sprintf("  Unique conds:   %d\n",
            length(unique(results_all$condition_id))))
cat(sprintf("  Successful fits: %d / %d\n",
            sum(results_all$status == "success", na.rm = TRUE),
            nrow(results_all)))
cat(sprintf("  Failed fits:     %d\n",
            sum(results_all$status == "failed", na.rm = TRUE)))
cat(sprintf("  Total time:      %s\n", format_time(overall_elapsed)))


# =============================================================================
# 6. Clean Up Parallel Backend
# =============================================================================

future::plan("sequential")
log_msg("Parallel backend shut down.", config = config)

cat("\nSimulation phase complete. Proceed to 03_run_postprocessing.R\n")
