# =============================================================================
# 02_run_simulation.R
# =============================================================================
# Phase 2: Main simulation loop.
#
# Strategy:
#   - Process CONDITIONS sequentially (memory management: one condition's
#     calibration result + datasets in memory at a time).
#   - Within each condition, REPLICATIONS are parallelized using
#     future + furrr (the unit of parallelization is run_one_replication).
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
# 1. Load Calibration Results
# =============================================================================

log_msg("Loading calibration results...", config = config)

calib_results <- vector("list", nrow(design))
names(calib_results) <- as.character(design$condition_id)

for (i in seq_len(nrow(design))) {
  cond_id <- design$condition_id[i]
  cache_path <- file.path(
    config$paths$calibration_dir,
    sprintf("calib_%04d.rds", cond_id)
  )
  if (!file.exists(cache_path)) {
    stop("Missing calibration for condition ", cond_id,
         ". Run 01_run_calibration.R first.")
  }
  calib_results[[i]] <- readRDS(cache_path)
}

log_msg("All ", nrow(design), " calibration results loaded.", config = config)


# =============================================================================
# 2. Configure Parallel Backend
# =============================================================================

n_workers <- config$parallel$n_workers
log_msg("Setting up parallel backend: ", config$parallel$strategy,
        " with ", n_workers, " workers.", config = config)

future::plan(config$parallel$strategy, workers = n_workers)

# Enable progressr globally
progressr::handlers("cli")


# =============================================================================
# 3. Main Simulation Loop (Conditions x Replications)
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

  design_row  <- design[ci, ]
  cond_id     <- design_row$condition_id
  calib_result <- calib_results[[ci]]

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
  progressr::with_progress({

    p <- progressr::progressor(along = remaining_reps)

    rep_results <- furrr::future_map(
      remaining_reps,
      function(rep_id) {

        # Run one complete replication
        result <- run_one_replication(
          design_row   = design_row,
          calib_result = calib_result,
          rep_id       = rep_id,
          config       = config
        )

        # Save checkpoint
        save_checkpoint(result, cond_id, rep_id, config)

        p(sprintf("cond %d rep %d", cond_id, rep_id))

        result
      },
      .options = furrr::furrr_options(
        seed = TRUE,              # Ensure RNG safety across workers
        chunk_size = config$parallel$chunk_size
      )
    )
  })

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
# 4. Compile All Results
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
# 5. Clean Up Parallel Backend
# =============================================================================

future::plan("sequential")
log_msg("Parallel backend shut down.", config = config)

cat("\nSimulation phase complete. Proceed to 03_run_postprocessing.R\n")
