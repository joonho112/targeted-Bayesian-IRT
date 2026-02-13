# =============================================================================
# 01_calibrate.R
# =============================================================================
# Pre-calibration phase: EQC calibration for each unique calibration cell.
#
# KEY DESIGN POINT: Calibration depends ONLY on (irt_model, latent_shape,
# target_rho) — NOT on N. There are only 30 unique calibration cells
# (2 models × 3 shapes × 5 rhos). The same calibrated item parameters are
# reused across all N levels and all K replications within a cell.
#
# Dependencies: IRTsimrel (>= 0.2.0)
#
# Author: JoonHo Lee
# Date: February 2026
# =============================================================================


#' Create Calibration Cell Grid
#'
#' Identifies the 30 unique (irt_model, latent_shape, target_rho)
#' combinations that require separate EQC calibrations.
#'
#' @param config Configuration list from load_config().
#' @return A tibble with columns: calib_cell_id, irt_model, latent_shape,
#'   I, target_rho.
#' @export
create_calibration_grid <- function(config) {

  grid <- tidyr::expand_grid(
    irt_model    = config$design$irt_models,
    latent_shape = config$design$latent_shapes,
    target_rho   = config$design$target_rhos
  )

  grid$I <- config$design$n_items

  # Create a human-readable cell ID
  grid$calib_cell_id <- sprintf(
    "%s_%s_rho%02d",
    grid$irt_model, grid$latent_shape,
    as.integer(grid$target_rho * 100)
  )

  grid <- dplyr::select(grid, calib_cell_id, irt_model, latent_shape, I, target_rho)
  grid
}


#' Run EQC Calibration for One Cell
#'
#' Performs Empirical Quadrature Calibration for a single
#' (irt_model, latent_shape, target_rho) combination.
#'
#' @param cell_row A one-row tibble from create_calibration_grid().
#' @param config Configuration list from load_config().
#' @return An eqc_result object from IRTsimrel::eqc_calibrate().
#' @export
run_calibration_cell <- function(cell_row, config) {

  cell_id <- cell_row$calib_cell_id

  # --- Check for cached result ---
  cache_path <- file.path(
    config$paths$calibration_dir,
    sprintf("eqc_%s.rds", cell_id)
  )

  if (file.exists(cache_path)) {
    log_msg("  [cache hit] ", cell_id, config = config)
    return(readRDS(cache_path))
  }

  # --- Derive deterministic seed for this calibration cell ---
  calib_seed <- derive_seed(
    master_seed  = config$seeds$master_seed,
    condition_id = paste0("calib_", cell_id)
  )

  # --- Run EQC calibration ---
  eqc_result <- IRTsimrel::eqc_calibrate(
    target_rho         = cell_row$target_rho,
    n_items            = as.integer(cell_row$I),
    model              = cell_row$irt_model,
    latent_shape       = cell_row$latent_shape,
    item_source        = config$calibration$item_source,
    reliability_metric = config$calibration$reliability_metric,
    M                  = as.integer(config$calibration$M_quadrature),
    c_bounds           = as.numeric(config$calibration$c_bounds),
    tol                = config$calibration$tol,
    seed               = calib_seed,
    verbose            = FALSE
  )

  # --- Verify achieved reliability ---
  rho_error <- abs(eqc_result$achieved_rho - cell_row$target_rho)

  log_msg(sprintf("  [%s] c*=%.4f, achieved_rho=%.4f (error=%.6f, status=%s)",
                  cell_id,
                  eqc_result$c_star,
                  eqc_result$achieved_rho,
                  rho_error,
                  eqc_result$misc$root_status),
          config = config)

  # --- Save to disk ---
  saveRDS(eqc_result, cache_path)

  eqc_result
}


#' Run Calibration for All 30 Cells
#'
#' @param config Configuration list from load_config().
#' @return A named list of eqc_result objects, keyed by calib_cell_id.
#' @export
run_all_calibrations <- function(config) {

  ensure_directories(config)

  grid <- create_calibration_grid(config)
  n_cells <- nrow(grid)

  log_msg(sprintf("=== Phase 1a: EQC Calibration (%d cells) ===", n_cells),
          config = config)

  calib_results <- vector("list", n_cells)
  names(calib_results) <- grid$calib_cell_id

  for (i in seq_len(n_cells)) {
    calib_results[[i]] <- run_calibration_cell(grid[i, ], config)
  }

  log_msg(sprintf("All %d calibration cells complete.", n_cells),
          config = config)

  calib_results
}


#' Look Up Calibration Result for a Design Condition
#'
#' Maps a condition (which includes N) to the correct calibration cell
#' (which does NOT depend on N).
#'
#' @param design_row A one-row tibble from create_design_matrix().
#' @param calib_results Named list from run_all_calibrations().
#' @return The eqc_result object for this condition's calibration cell.
#' @export
lookup_calibration <- function(design_row, calib_results) {

  cell_id <- sprintf(
    "%s_%s_rho%02d",
    design_row$irt_model,
    design_row$latent_shape,
    as.integer(design_row$target_rho * 100)
  )

  if (!cell_id %in% names(calib_results)) {
    stop("Calibration cell not found: ", cell_id, call. = FALSE)
  }

  calib_results[[cell_id]]
}


#' Build Calibration Summary Table
#'
#' Creates a tidy data frame summarizing all calibration results.
#'
#' @param calib_results Named list from run_all_calibrations().
#' @return A tibble with one row per calibration cell.
#' @export
build_calibration_summary <- function(calib_results) {

  rows <- lapply(names(calib_results), function(cell_id) {
    res <- calib_results[[cell_id]]
    data.frame(
      calib_cell_id = cell_id,
      irt_model     = res$model,
      n_items       = res$n_items,
      target_rho    = res$target_rho,
      achieved_rho  = res$achieved_rho,
      c_star        = res$c_star,
      rho_error     = abs(res$achieved_rho - res$target_rho),
      root_status   = res$misc$root_status,
      beta_mean     = mean(res$beta_vec),
      beta_sd       = sd(res$beta_vec),
      lambda_mean   = mean(res$lambda_scaled),
      lambda_sd     = sd(res$lambda_scaled),
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(rows)
}
