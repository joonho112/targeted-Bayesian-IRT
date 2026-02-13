# =============================================================================
# 04_evaluate.R
# =============================================================================
# Replication-level orchestration and result assembly.
#
# The main entry point is run_one_replication() which:
#   1. Generates one dataset
#   2. Fits all 3 priors
#   3. Collects losses into a tidy data frame
#   4. Saves checkpoint
#
# Also contains: standalone loss functions, checkpoint utilities,
# result compilation.
#
# Author: JoonHo Lee
# Date: February 2026
# =============================================================================


# =============================================================================
# Standalone Loss Functions (for validation)
# =============================================================================

#' @export
compute_rmse <- function(est, true) {
  sqrt(mean((est - true)^2))
}

#' @export
compute_mselp <- function(est, true) {
  N <- length(est)
  mean((rank(est) / N - rank(true) / N)^2)
}

#' @export
compute_ks <- function(est, true) {
  suppressWarnings(ks.test(est, true)$statistic)
}


# =============================================================================
# Replication-Level Orchestration
# =============================================================================

#' Run One Complete Replication
#'
#' The primary unit of work. For one (condition, replication):
#'   1. Generate dataset via simulate_response_data()
#'   2. Fit all 3 priors via dpmirt()
#'   3. Compute PM/CB/GR estimates and losses
#'   4. Return tidy results tibble + save checkpoint
#'
#' @param design_row One-row tibble from design matrix.
#' @param eqc_result Calibration result for this condition's cell.
#' @param rep_id Integer. Replication index (1..K).
#' @param alpha_cache Named list from build_alpha_cache().
#' @param config Configuration list.
#'
#' @return Tibble with one row per (prior x summary_method x parameter).
#' @export
run_one_replication <- function(design_row, eqc_result, rep_id,
                                alpha_cache, config) {

  cond_id <- design_row$condition_id

  # --- Step 1: Generate dataset ---
  rep_seed <- derive_seed(config$seeds$master_seed, cond_id, rep_id = rep_id)
  dataset <- generate_one_dataset(
    eqc_result   = eqc_result,
    N            = design_row$N,
    latent_shape = design_row$latent_shape,
    rep_seed     = rep_seed
  )

  # --- Step 2: Fit all 3 priors ---
  fit_results <- fit_all_priors(
    dataset        = dataset,
    design_row     = design_row,
    rep_id         = rep_id,
    alpha_cache    = alpha_cache,
    config         = config,
    save_posteriors = config$posterior_storage$save_posteriors
  )

  # --- Step 3: Assemble tidy results ---
  rows <- list()
  for (prior_name in names(fit_results$prior_results)) {
    pr <- fit_results$prior_results[[prior_name]]

    if (is.null(pr$losses)) {
      # Failed fit - create placeholder
      rows[[prior_name]] <- data.frame(
        condition_id = cond_id,
        irt_model    = design_row$irt_model,
        latent_shape = design_row$latent_shape,
        N            = design_row$N,
        I            = design_row$I,
        target_rho   = design_row$target_rho,
        rep_id       = rep_id,
        prior        = prior_name,
        parameter    = "theta",
        method       = "pm",
        msel         = NA_real_,
        mselr        = NA_real_,
        ks           = NA_real_,
        status       = "failed",
        time_secs    = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }

    loss_df <- pr$losses
    loss_df$condition_id <- cond_id
    loss_df$irt_model    <- design_row$irt_model
    loss_df$latent_shape <- design_row$latent_shape
    loss_df$N            <- design_row$N
    loss_df$I            <- design_row$I
    loss_df$target_rho   <- design_row$target_rho
    loss_df$rep_id       <- rep_id
    loss_df$prior        <- prior_name
    loss_df$status       <- pr$fit_summary$status
    loss_df$time_secs    <- pr$fit_summary$time_secs

    rows[[prior_name]] <- loss_df
  }

  result_df <- dplyr::bind_rows(rows)

  # --- Step 4: Save checkpoint ---
  if (isTRUE(config$checkpointing$enabled)) {
    save_checkpoint(result_df, cond_id, rep_id, config)
  }

  result_df
}


# =============================================================================
# Checkpoint Utilities
# =============================================================================

#' @export
save_checkpoint <- function(result_df, cond_id, rep_id, config) {
  dir.create(config$paths$checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(
    config$paths$checkpoint_dir,
    sprintf("ckpt_c%04d_r%03d.rds", cond_id, rep_id)
  )
  saveRDS(result_df, path)
  invisible(path)
}

#' @export
get_completed_reps <- function(cond_id, config) {
  ckpt_dir <- config$paths$checkpoint_dir
  if (!dir.exists(ckpt_dir)) return(integer(0))

  pattern <- sprintf("ckpt_c%04d_r(\\d+)\\.rds", cond_id)
  files <- list.files(ckpt_dir, pattern = pattern)
  if (length(files) == 0) return(integer(0))

  m <- regmatches(files, regexec("_r(\\d+)\\.rds", files))
  sort(as.integer(vapply(m, function(x) x[2], character(1))))
}

#' @export
compile_condition_results <- function(cond_id, config) {
  pattern <- sprintf("ckpt_c%04d_r\\d+\\.rds", cond_id)
  files <- list.files(config$paths$checkpoint_dir, pattern = pattern,
                      full.names = TRUE)
  if (length(files) == 0) return(tibble::tibble())
  dplyr::bind_rows(lapply(files, readRDS))
}

#' @export
compile_all_results <- function(config) {
  files <- list.files(config$paths$checkpoint_dir, pattern = "ckpt_.*\\.rds",
                      full.names = TRUE)
  if (length(files) == 0) stop("No checkpoint files found.")
  dplyr::bind_rows(lapply(files, readRDS))
}
