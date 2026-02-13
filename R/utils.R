# =============================================================================
# utils.R
# =============================================================================
# Shared helper functions used across the simulation pipeline.
#
# Contents:
#   - safe_require(): Safely check for package availability
#   - format_time(): Human-readable time formatting
#   - format_bytes(): Human-readable file size formatting
#   - summarize_results(): Quick summary statistics of simulation results
#   - estimate_storage(): Estimate total posterior storage requirements
#   - validate_design(): Validate design matrix integrity
#
# Author: JoonHo Lee
# Date: February 2026
# =============================================================================


#' Safely Require a Package
#'
#' Checks if a package is installed and loadable. Provides a helpful
#' error message if not.
#'
#' @param pkg Character. Package name.
#' @param purpose Character. Brief description of why the package is needed.
#'
#' @return TRUE if package is available.
#' @export
safe_require <- function(pkg, purpose = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    msg <- sprintf("Package '%s' is required", pkg)
    if (!is.null(purpose)) {
      msg <- paste0(msg, " for ", purpose)
    }
    msg <- paste0(msg, ". Install with: install.packages('", pkg, "')")
    stop(msg, call. = FALSE)
  }
  invisible(TRUE)
}


#' Format Elapsed Time as Human-Readable String
#'
#' @param seconds Numeric. Elapsed time in seconds.
#' @return Character string (e.g., "2h 15m 30s").
#' @export
format_time <- function(seconds) {
  if (is.na(seconds)) return("NA")

  hours   <- floor(seconds / 3600)
  minutes <- floor((seconds %% 3600) / 60)
  secs    <- round(seconds %% 60, 1)

  if (hours > 0) {
    sprintf("%dh %dm %gs", hours, minutes, secs)
  } else if (minutes > 0) {
    sprintf("%dm %gs", minutes, secs)
  } else {
    sprintf("%.1fs", secs)
  }
}


#' Format Bytes as Human-Readable String
#'
#' @param bytes Numeric. Size in bytes.
#' @return Character string (e.g., "1.23 GB").
#' @export
format_bytes <- function(bytes) {
  if (bytes < 1024) {
    sprintf("%d B", bytes)
  } else if (bytes < 1024^2) {
    sprintf("%.1f KB", bytes / 1024)
  } else if (bytes < 1024^3) {
    sprintf("%.1f MB", bytes / 1024^2)
  } else {
    sprintf("%.2f GB", bytes / 1024^3)
  }
}


#' Estimate Total Posterior Storage Requirements
#'
#' Computes the expected disk usage for saving all posterior draws across
#' the full factorial design.
#'
#' @param config A configuration list from load_config().
#'
#' @return A list with per_file_bytes (by N), total_bytes, total_human.
#' @export
estimate_storage <- function(config) {

  n_retained <- config$mcmc$niter - config$mcmc$nburnin
  n_priors   <- 3  # Gaussian, DP-Broad, DP-Focused
  K          <- config$design$replications

  # Bytes per element (double precision)
  bytes_per_element <- 8

  # Number of conditions
  n_conditions <- length(config$design$irt_models) *
    length(config$design$latent_shapes) *
    length(config$design$sample_sizes) *
    length(config$design$target_rhos)

  # Per-file size depends on N (person count)
  per_N_sizes <- sapply(config$design$sample_sizes, function(N) {
    # theta_samp: n_retained x N
    # beta_samp:  n_retained x I
    # lambda_samp: n_retained x I (for 2PL; NULL for Rasch)
    I <- config$design$n_items

    theta_bytes  <- n_retained * N * bytes_per_element
    beta_bytes   <- n_retained * I * bytes_per_element
    lambda_bytes <- n_retained * I * bytes_per_element

    # Approximate: compressed RDS is ~30-50% of raw
    raw_bytes <- theta_bytes + beta_bytes + lambda_bytes
    compressed_bytes <- raw_bytes * 0.4  # Conservative estimate

    c(raw = raw_bytes, compressed = compressed_bytes)
  })

  colnames(per_N_sizes) <- config$design$sample_sizes

  # Total files: conditions * priors * replications
  # But conditions are spread across N values, so compute per-N
  n_cond_per_N <- n_conditions / length(config$design$sample_sizes)
  total_files <- n_conditions * n_priors * K

  total_bytes <- sum(sapply(config$design$sample_sizes, function(N) {
    n_cond_per_N * n_priors * K * per_N_sizes["compressed", as.character(N)]
  }))

  list(
    per_file_sizes = per_N_sizes,
    total_files    = total_files,
    total_bytes    = total_bytes,
    total_human    = format_bytes(total_bytes)
  )
}


#' Validate Design Matrix Integrity
#'
#' Checks that the design matrix is well-formed and consistent with config.
#'
#' @param design A tibble from create_design_matrix().
#' @param config A configuration list.
#'
#' @return TRUE if valid; errors otherwise.
#' @export
validate_design <- function(design, config) {

  # Check required columns
  required_cols <- c("condition_id", "irt_model", "latent_shape",
                     "N", "I", "target_rho")
  missing <- setdiff(required_cols, names(design))
  if (length(missing) > 0) {
    stop("Design matrix missing columns: ", paste(missing, collapse = ", "))
  }

  # Check no duplicate condition_ids
  if (any(duplicated(design$condition_id))) {
    stop("Duplicate condition_ids found in design matrix.")
  }

  # Check expected dimensions
  expected_rows <- length(config$design$irt_models) *
    length(config$design$latent_shapes) *
    length(config$design$sample_sizes) *
    length(config$design$target_rhos)

  if (nrow(design) != expected_rows) {
    warning(
      "Design has ", nrow(design), " rows; expected ", expected_rows,
      " from full factorial cross."
    )
  }

  # Check I is consistent
  if (!all(design$I == config$design$n_items)) {
    stop("Inconsistent I values in design matrix.")
  }

  invisible(TRUE)
}


#' Summarize Simulation Results
#'
#' Produces quick summary statistics (mean, median, SD) of loss metrics
#' across replications for each condition.
#'
#' @param results_df A tibble of compiled simulation results.
#' @param trim Numeric. Fraction to trim when computing trimmed mean.
#'   Default 0 (no trimming).
#'
#' @return A tibble of summary statistics per condition.
#' @export
summarize_results <- function(results_df, trim = 0) {

  metric_cols <- intersect(
    c("rmse", "mselp", "ks", "msel", "mselr"),
    names(results_df)
  )

  group_cols <- c("condition_id", "irt_model", "latent_shape",
                  "N", "I", "target_rho", "prior", "summary_method",
                  "parameter")

  # Only include existing columns in grouping
  group_cols <- intersect(group_cols, names(results_df))

  summary_df <- results_df |>
    dplyr::filter(status == "success") |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarize(
      n_reps = dplyr::n(),
      dplyr::across(
        dplyr::all_of(metric_cols),
        list(
          mean   = ~ mean(.x, na.rm = TRUE, trim = trim),
          median = ~ median(.x, na.rm = TRUE),
          sd     = ~ sd(.x, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    )

  summary_df
}


#' Print a Concise Progress Summary
#'
#' @param design A tibble from create_design_matrix().
#' @param config A configuration list.
#'
#' @return Invisible NULL. Prints to console.
#' @export
print_study_summary <- function(design, config) {

  n_conditions   <- nrow(design)
  n_priors       <- 3
  K              <- config$design$replications
  n_fits_total   <- n_conditions * n_priors * K
  n_workers      <- config$parallel$n_workers

  cat("\n")
  cat("=======================================================\n")
  cat("  Targeted Bayesian IRT Simulation Study\n")
  cat("=======================================================\n\n")
  cat(sprintf("  Conditions:       %d\n", n_conditions))
  cat(sprintf("  Priors per cond:  %d\n", n_priors))
  cat(sprintf("  Replications:     %d\n", K))
  cat(sprintf("  Total model fits: %d\n", n_fits_total))
  cat(sprintf("  Parallel workers: %d\n", n_workers))
  cat(sprintf("  Master seed:      %d\n", config$seeds$master_seed))
  cat("\n")

  storage <- estimate_storage(config)
  cat(sprintf("  Est. posterior storage: %s (%d files)\n",
              storage$total_human, storage$total_files))
  cat("\n")
  cat("  Design factors:\n")
  cat(sprintf("    IRT models:   %s\n",
              paste(config$design$irt_models, collapse = ", ")))
  cat(sprintf("    Latent G:     %s\n",
              paste(config$design$latent_shapes, collapse = ", ")))
  cat(sprintf("    N:            %s\n",
              paste(config$design$sample_sizes, collapse = ", ")))
  cat(sprintf("    I:            %d\n", config$design$n_items))
  cat(sprintf("    Target rho:   %s\n",
              paste(config$design$target_rhos, collapse = ", ")))
  cat("\n")
  cat("=======================================================\n\n")

  invisible(NULL)
}
