# =============================================================================
# 03_fit.R
# =============================================================================
# Model fitting module: fit Bayesian IRT models with DPMirt.
#
# For each dataset x prior combination:
#   1. Fit the model via DPMirt::dpmirt()
#   2. Compute PM/CB/GR estimates via DPMirt::dpmirt_estimates()
#   3. Compute loss metrics via DPMirt::dpmirt_loss()
#   4. Optionally save posterior draws to disk
#
# IMPORTANT: library(nimble) must be loaded before library(DPMirt)
#
# Dependencies: DPMirt (>= 0.1.0), DPprior (>= 0.1.0), nimble
#
# Author: JoonHo Lee
# Date: February 2026
# =============================================================================


#' Build Alpha Prior Cache
#'
#' Pre-computes the Gamma(a,b) hyperprior for each DP variant x N combination.
#' Only 8 unique calls: 2 DP priors x 4 N levels.
#'
#' @param config Configuration list.
#' @return Named list of c(a, b) vectors, keyed by "{prior_key}_N{N}".
#' @export
build_alpha_cache <- function(config) {

  cache <- list()

  dp_priors <- list(
    dp_broad   = config$priors$dp_broad,
    dp_focused = config$priors$dp_focused
  )

  for (prior_name in names(dp_priors)) {
    prior_cfg <- dp_priors[[prior_name]]

    for (N in config$design$sample_sizes) {
      mu_K <- max(prior_cfg$mu_K_fraction * N, 2)  # mu_K must be > 1

      dp_fit <- DPprior::DPprior_fit(
        J          = as.integer(N),
        mu_K       = mu_K,
        confidence = prior_cfg$confidence,
        method     = "A2-MN",
        check_diagnostics = FALSE,
        warn_dominance    = FALSE,
        verbose           = FALSE
      )

      cache_key <- sprintf("%s_N%d", prior_name, N)
      cache[[cache_key]] <- c(a = dp_fit$a, b = dp_fit$b)
    }
  }

  cache
}


#' Fit One Model to One Dataset
#'
#' @param dataset List from generate_one_dataset() with response_matrix, theta, beta, lambda.
#' @param irt_model Character. "rasch" or "2pl".
#' @param prior_spec List with: label, type ("normal"/"dpm"), prior_key, prior_index.
#'   For DPM: also mu_K, confidence.
#' @param alpha_prior Numeric c(a,b) or NULL. Pre-computed Gamma hyperprior.
#' @param model_seed Integer. MCMC seed.
#' @param config Configuration list.
#'
#' @return List with: fit, estimates, losses, fit_summary, prior_key.
#' @export
fit_one_model <- function(dataset, irt_model, prior_spec, alpha_prior,
                          model_seed, config) {

  t_start <- Sys.time()

  # --- Build dpmirt() arguments ---
  dpmirt_args <- list(
    data       = dataset$response_matrix,
    model      = irt_model,
    prior      = prior_spec$type,
    niter      = as.integer(config$mcmc$niter),
    nburnin    = as.integer(config$mcmc$nburnin),
    thin       = as.integer(config$mcmc$thin),
    nchains    = as.integer(config$mcmc$nchains),
    seed       = model_seed,
    rescale    = TRUE,
    compute_waic       = FALSE,
    compute_dp_density = FALSE,
    save_draws = TRUE,
    verbose    = FALSE
  )

  # Add DPM-specific settings
  if (prior_spec$type == "dpm") {
    dpmirt_args$alpha_prior <- alpha_prior
    dpmirt_args$M           <- as.integer(config$mcmc$M_dpm)
  }

  # --- Fit model ---
  fit <- tryCatch(
    do.call(DPMirt::dpmirt, dpmirt_args),
    error = function(e) {
      warning(sprintf("  MCMC FAILED [%s]: %s", prior_spec$label, conditionMessage(e)))
      return(NULL)
    }
  )

  if (is.null(fit)) {
    return(list(
      fit         = NULL,
      estimates   = NULL,
      losses      = NULL,
      fit_summary = list(status = "failed", error = "dpmirt() error"),
      prior_key   = prior_spec$prior_key
    ))
  }

  # --- Compute PM/CB/GR estimates ---
  estimates <- DPMirt::dpmirt_estimates(
    fit,
    methods = c("pm", "cb", "gr"),
    alpha   = 0.05
  )

  # --- Compute loss metrics ---
  true_lambda <- if (irt_model == "2pl") dataset$lambda else NULL

  losses <- DPMirt::dpmirt_loss(
    estimates  = estimates,
    true_theta = dataset$theta,
    true_beta  = dataset$beta,
    true_lambda = true_lambda,
    metrics    = c("msel", "mselr", "ks")
  )

  t_end <- Sys.time()
  elapsed <- as.numeric(difftime(t_end, t_start, units = "secs"))

  # --- Build fit summary ---
  fit_summary <- list(
    status    = "success",
    time_secs = elapsed,
    n_post    = nrow(fit$theta_samp),
    N         = ncol(fit$theta_samp)
  )

  list(
    fit         = fit,
    estimates   = estimates,
    losses      = losses,
    fit_summary = fit_summary,
    prior_key   = prior_spec$prior_key
  )
}


#' Fit All Three Priors for One Dataset
#'
#' @param dataset List from generate_one_dataset().
#' @param design_row One-row tibble from design matrix.
#' @param rep_id Integer. Replication index.
#' @param alpha_cache Named list from build_alpha_cache().
#' @param config Configuration list.
#' @param save_posteriors Logical. Whether to save posterior draws.
#'
#' @return List of 3 fit results (one per prior), plus metadata.
#' @export
fit_all_priors <- function(dataset, design_row, rep_id,
                           alpha_cache, config,
                           save_posteriors = config$posterior_storage$save_posteriors) {

  prior_specs <- get_prior_specs(config, design_row$N)
  master_seed <- config$seeds$master_seed
  cond_id     <- design_row$condition_id
  irt_model   <- design_row$irt_model

  results <- vector("list", length(prior_specs))
  names(results) <- names(prior_specs)

  for (prior_name in names(prior_specs)) {
    pspec <- prior_specs[[prior_name]]

    # Derive model seed
    model_seed <- derive_seed(master_seed, cond_id,
                              rep_id = rep_id,
                              prior_index = pspec$prior_index)

    # Look up alpha prior from cache (NULL for Gaussian)
    alpha_prior <- NULL
    if (pspec$type == "dpm") {
      cache_key <- sprintf("%s_N%d", prior_name, design_row$N)
      alpha_prior <- alpha_cache[[cache_key]]
    }

    # Fit
    result <- fit_one_model(
      dataset     = dataset,
      irt_model   = irt_model,
      prior_spec  = pspec,
      alpha_prior = alpha_prior,
      model_seed  = model_seed,
      config      = config
    )

    # Optionally save posteriors
    if (save_posteriors && !is.null(result$fit)) {
      posterior_path <- save_posterior_draws(
        fit       = result$fit,
        cond_id   = cond_id,
        prior_key = prior_name,
        rep_id    = rep_id,
        config    = config
      )
      result$posteriors_path <- posterior_path
    }

    # Free MCMC fit object to save memory (keep estimates and losses)
    result$fit <- NULL

    results[[prior_name]] <- result
  }

  list(
    condition_id = cond_id,
    rep_id       = rep_id,
    irt_model    = irt_model,
    N            = design_row$N,
    prior_results = results
  )
}


#' Save Posterior Draws to Disk
#'
#' @param fit dpmirt_fit object.
#' @param cond_id Integer or character. Condition ID.
#' @param prior_key Character. Prior key name.
#' @param rep_id Integer. Replication ID.
#' @param config Configuration list.
#'
#' @return File path (character).
#' @export
save_posterior_draws <- function(fit, cond_id, prior_key, rep_id, config) {

  dir_path <- file.path(
    config$paths$posteriors_dir,
    sprintf("cond_%04d", as.integer(cond_id)),
    prior_key
  )
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)

  file_path <- file.path(dir_path, sprintf("rep_%03d.rds", rep_id))

  posterior_data <- list(
    theta_samp  = fit$theta_samp,
    beta_samp   = fit$beta_samp,
    lambda_samp = fit$lambda_samp
  )

  saveRDS(posterior_data, file_path, compress = TRUE)
  file_path
}
