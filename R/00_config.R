# =============================================================================
# 00_config.R
# =============================================================================
# Configuration loader, design matrix constructor, and seed management.
#
# Contents:
#   - load_config(): Load and validate YAML configuration
#   - create_design_matrix(): Full factorial design as a tibble
#   - get_prior_specs(): Build prior specification list for a given N
#   - derive_seed(): Deterministic seed derivation
#   - ensure_directories(): Create output directory tree
#
# Dependencies: yaml, tibble, dplyr, digest
#
# Author: JoonHo Lee
# Date: February 2026
# =============================================================================


#' Load and Validate YAML Configuration
#'
#' Reads the master configuration file and resolves paths relative to the
#' project root.
#'
#' @param config_path Character. Path to YAML configuration file.
#'   Default: "config/sim_config.yaml" relative to project root.
#' @param project_root Character. Project root directory. If NULL, inferred
#'   from the location of this source file or the current working directory.
#'
#' @return A validated configuration list.
#'
#' @export
load_config <- function(config_path = NULL, project_root = NULL) {

  # --- Infer project root ---
  if (is.null(project_root)) {
    # Try rstudioapi, then fall back to getwd()
    project_root <- tryCatch(
      rstudioapi::getActiveProject(),
      error = function(e) getwd()
    )
    if (is.null(project_root)) project_root <- getwd()
  }

  # --- Resolve config path ---
  if (is.null(config_path)) {
    config_path <- file.path(project_root, "config", "sim_config.yaml")
  }

  if (!file.exists(config_path)) {
    stop("Configuration file not found: ", config_path, call. = FALSE)
  }

  config <- yaml::read_yaml(config_path)

  # --- Resolve paths ---
  config$paths$project_root <- normalizePath(project_root, mustWork = TRUE)

  path_fields <- c("output_dir", "calibration_dir", "posteriors_dir",
                    "results_dir", "figures_dir", "checkpoint_dir")
  for (field in path_fields) {
    config$paths[[field]] <- file.path(
      config$paths$project_root,
      config$paths[[field]]
    )
  }

  # --- Validate critical settings ---
  stopifnot(
    is.numeric(config$design$n_items),
    config$design$n_items > 0,
    is.numeric(config$design$replications),
    config$design$replications > 0,
    all(config$design$target_rhos > 0 & config$design$target_rhos < 1),
    all(config$design$sample_sizes > 0)
  )

  config
}


#' Create Full Factorial Design Matrix
#'
#' Generates a tibble where each row represents one unique experimental
#' condition. The design crosses: irt_model x latent_shape x N x target_rho.
#'
#' @param config A configuration list from load_config().
#'
#' @return A tibble with columns:
#'   condition_id (integer), irt_model (character), latent_shape (character),
#'   N (integer), I (integer), target_rho (numeric).
#'
#' @export
create_design_matrix <- function(config) {

  design <- tidyr::expand_grid(
    irt_model    = config$design$irt_models,
    latent_shape = config$design$latent_shapes,
    N            = config$design$sample_sizes,
    target_rho   = config$design$target_rhos
  )

  design <- dplyr::mutate(
    design,
    I            = config$design$n_items,
    condition_id = dplyr::row_number()
  )

  # Reorder columns for clarity
  design <- dplyr::select(
    design,
    condition_id, irt_model, latent_shape, N, I, target_rho
  )

  design
}


#' Build Prior Specification List for a Given Sample Size
#'
#' Translates the YAML prior configuration into a list of prior
#' specification lists, resolving the mu_K_fraction * N computation.
#'
#' @param config A configuration list from load_config().
#' @param N Integer. Current sample size (needed for DP mu_K computation).
#'
#' @return A named list of prior specifications. Each element is a list
#'   with fields: label, type, mu_K (for DPM), confidence (for DPM).
#'
#' @export
get_prior_specs <- function(config, N) {

  prior_specs <- list()

  # --- Gaussian prior ---
  prior_specs[["gaussian"]] <- list(
    label      = config$priors$gaussian$label,
    type       = config$priors$gaussian$type,
    prior_key  = "gaussian",
    prior_index = 1L
  )

  # --- DP-Broad prior ---
  dp_broad_mu_K <- config$priors$dp_broad$mu_K_fraction * N
  # mu_K must be > 1 for DPprior_fit

  dp_broad_mu_K <- max(dp_broad_mu_K, 2)
  prior_specs[["dp_broad"]] <- list(
    label       = config$priors$dp_broad$label,
    type        = config$priors$dp_broad$type,
    mu_K        = dp_broad_mu_K,
    confidence  = config$priors$dp_broad$confidence,
    prior_key   = "dp_broad",
    prior_index = 2L
  )

  # --- DP-Focused prior ---
  dp_focused_mu_K <- config$priors$dp_focused$mu_K_fraction * N
  dp_focused_mu_K <- max(dp_focused_mu_K, 2)
  prior_specs[["dp_focused"]] <- list(
    label       = config$priors$dp_focused$label,
    type        = config$priors$dp_focused$type,
    mu_K        = dp_focused_mu_K,
    confidence  = config$priors$dp_focused$confidence,
    prior_key   = "dp_focused",
    prior_index = 3L
  )

  prior_specs
}


# =============================================================================
# Seed Management
# =============================================================================

#' Derive a Deterministic Seed
#'
#' Produces reproducible, unique seeds for each level of the simulation
#' hierarchy using hashing. The hierarchy is:
#'
#'   master_seed -> condition_seed -> replication_seed -> model_seed
#'
#' @param master_seed Integer. The global master seed from config.
#' @param condition_id Integer. Unique condition identifier.
#' @param rep_id Integer or NULL. Replication identifier.
#' @param prior_index Integer or NULL. Prior model index (1, 2, or 3).
#'
#' @return Integer. A derived seed value.
#'
#' @details
#' Uses digest::digest2int for deterministic hashing. The seed is always
#' positive and within .Machine$integer.max.
#'
#' @export
derive_seed <- function(master_seed, condition_id,
                        rep_id = NULL, prior_index = NULL) {

  # Condition-level seed
  seed_string <- paste(master_seed, condition_id, sep = "_")

  if (!is.null(rep_id)) {
    seed_string <- paste(seed_string, rep_id, sep = "_")
  }

  if (!is.null(prior_index)) {
    seed_string <- paste(seed_string, prior_index, sep = "_")
  }

  # Use digest to produce a deterministic integer
  raw_hash <- digest::digest2int(seed_string)

  # Ensure positive seed within valid range
  abs(raw_hash) %% (.Machine$integer.max - 1L) + 1L
}


# =============================================================================
# Directory Management
# =============================================================================

#' Ensure All Output Directories Exist
#'
#' Creates the full output directory tree as specified in config.
#'
#' @param config A configuration list from load_config().
#'
#' @return Invisible NULL. Called for side effect of creating directories.
#'
#' @export
ensure_directories <- function(config) {

  dirs <- c(
    config$paths$calibration_dir,
    config$paths$posteriors_dir,
    config$paths$results_dir,
    config$paths$figures_dir,
    config$paths$checkpoint_dir
  )

  for (d in dirs) {
    if (!dir.exists(d)) {
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
    }
  }

  invisible(NULL)
}


# =============================================================================
# Logging Utility
# =============================================================================

#' Log a Message with Timestamp
#'
#' Prints to console and optionally writes to log file.
#'
#' @param ... Arguments passed to paste0() to form the message.
#' @param config Optional configuration list (for log file path).
#' @param level Character. Log level: "INFO", "WARN", "ERROR".
#'
#' @return Invisible NULL.
#'
#' @export
log_msg <- function(..., config = NULL, level = "INFO") {

  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  msg <- paste0("[", timestamp, "] [", level, "] ", paste0(...))

  message(msg)

  # Write to log file if config provided
  if (!is.null(config) && !is.null(config$logging$log_file)) {
    log_path <- file.path(config$paths$project_root,
                          config$logging$log_file)
    cat(msg, "\n", file = log_path, append = TRUE)
  }

  invisible(NULL)
}
