# =============================================================================
# 02_generate.R
# =============================================================================
# Data generation module: generates response data for each (condition, rep).
#
# Uses IRTsimrel::simulate_response_data() which takes the pre-calibrated
# eqc_result and generates theta + response matrix consistently.
#
# KEY: Item parameters (beta, lambda) are FIXED from calibration.
#      Only theta and Y vary across replications.
#
# Dependencies: IRTsimrel (>= 0.2.0)
#
# Author: JoonHo Lee
# Date: February 2026
# =============================================================================


#' Generate One Dataset for a Single Replication
#'
#' Uses IRTsimrel::simulate_response_data() to draw fresh person abilities
#' and generate responses using the calibrated item parameters.
#'
#' @param eqc_result An eqc_result object from run_calibration_cell().
#' @param N Integer. Number of examinees.
#' @param latent_shape Character. Latent distribution shape.
#' @param rep_seed Integer. Deterministic seed for this replication.
#'
#' @return A list with components:
#'   \item{response_matrix}{N x I binary response matrix}
#'   \item{theta}{Numeric vector of true person abilities (length N)}
#'   \item{beta}{Numeric vector of item difficulties (length I)}
#'   \item{lambda}{Numeric vector of item discriminations (length I)}
#'   \item{N}{Sample size}
#'   \item{I}{Number of items}
#'
#' @export
generate_one_dataset <- function(eqc_result, N, latent_shape, rep_seed) {

  sim_data <- IRTsimrel::simulate_response_data(
    result       = eqc_result,
    n_persons    = as.integer(N),
    latent_shape = latent_shape,
    seed         = rep_seed
  )

  list(
    response_matrix = sim_data$response_matrix,
    theta           = sim_data$theta,
    beta            = sim_data$beta,
    lambda          = sim_data$lambda,
    N               = as.integer(N),
    I               = as.integer(ncol(sim_data$response_matrix))
  )
}


#' Generate K Datasets for One Condition
#'
#' @param design_row One-row tibble from create_design_matrix().
#' @param eqc_result EQC calibration result for this condition's cell.
#' @param config Configuration list.
#' @param K Number of replications (default: from config).
#'
#' @return List of K dataset lists.
#' @export
generate_condition_datasets <- function(design_row, eqc_result, config,
                                        K = config$design$replications) {

  master_seed <- config$seeds$master_seed
  cond_id     <- design_row$condition_id

  datasets <- vector("list", K)

  for (k in seq_len(K)) {
    rep_seed <- derive_seed(master_seed, cond_id, rep_id = k)
    datasets[[k]] <- generate_one_dataset(
      eqc_result   = eqc_result,
      N            = design_row$N,
      latent_shape = design_row$latent_shape,
      rep_seed     = rep_seed
    )
  }

  datasets
}
