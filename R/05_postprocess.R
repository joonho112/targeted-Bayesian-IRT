# =============================================================================
# 05_postprocess.R
# =============================================================================
# Post-processing: Scaled EDF, meta-model regression, and visualization.
#
# Contents:
#   - compute_scaled_edf(): Pooled binned density across replications
#   - run_meta_regression(): Log-linear regression with cluster-robust SEs
#   - Various plotting functions
#
# Dependencies: dplyr, tidyr, ggplot2, sandwich, lmtest, IRTsimrel
#
# Author: JoonHo Lee
# Date: February 2026
# =============================================================================


# =============================================================================
# Scaled Empirical Distribution Function
# =============================================================================

#' Compute Scaled EDF from Pooled Point Estimates
#'
#' Pools point estimates (PM, CB, or GR) across K replications and computes
#' binned density. Normalization: density = count / (K * N * bin_width).
#' This matches the old codebase pattern.
#'
#' @param est_list List of K numeric vectors (one per replication).
#' @param true_list List of K numeric vectors (true theta per rep).
#' @param latent_shape Character. For super-population reference.
#' @param n_bins Integer. Number of bins (default: 50).
#' @param theta_range Numeric length-2. Bin range (default: c(-6, 6)).
#' @param n_superpop Integer. Draws for super-population density.
#'
#' @return Data frame with columns: bin_mid, dens_est, dens_true_fs, dens_true_sp.
#' @export
compute_scaled_edf <- function(est_list, true_list, latent_shape,
                                n_bins = 50L, theta_range = c(-6, 6),
                                n_superpop = 100000L) {

  K <- length(est_list)
  N <- length(est_list[[1]])

  # Define bins
  breaks    <- seq(theta_range[1], theta_range[2], length.out = n_bins + 1)
  bin_width <- diff(breaks)[1]
  bin_mids  <- (breaks[-length(breaks)] + breaks[-1]) / 2

  # Pool estimates and true theta across replications
  pooled_est  <- unlist(est_list)
  pooled_true <- unlist(true_list)

  # Compute binned counts
  counts_est  <- as.numeric(table(cut(pooled_est,  breaks, include.lowest = TRUE)))
  counts_true <- as.numeric(table(cut(pooled_true, breaks, include.lowest = TRUE)))

  # Normalize: density = count / (total * bin_width)
  total <- K * N
  dens_est     <- counts_est  / (total * bin_width)
  dens_true_fs <- counts_true / (total * bin_width)

  # Super-population reference density
  sp_draws <- IRTsimrel::sim_latentG(
    n = n_superpop, shape = latent_shape, seed = 12345L
  )$theta
  counts_sp <- as.numeric(table(cut(sp_draws, breaks, include.lowest = TRUE)))
  dens_true_sp <- counts_sp / (n_superpop * bin_width)

  data.frame(
    bin_mid      = bin_mids,
    dens_est     = dens_est,
    dens_true_fs = dens_true_fs,
    dens_true_sp = dens_true_sp
  )
}


# =============================================================================
# Meta-Model Regression
# =============================================================================

#' Run Meta-Model Regression
#'
#' Fits: log(metric) ~ (N + rho + prior + method)^2
#' with cluster-robust standard errors (clustered on condition_id).
#'
#' @param results_df Tibble of compiled results (theta rows only).
#' @param metric Character. Column name of the metric to model.
#' @param config Configuration list.
#'
#' @return List with: model, robust_coefs, r_squared, formula_used.
#' @export
run_meta_regression <- function(results_df, metric = "msel", config = NULL) {

  # Prepare data
  df <- results_df %>%
    dplyr::filter(parameter == "theta", !is.na(.data[[metric]])) %>%
    dplyr::mutate(
      log_metric = log(.data[[metric]]),
      N_fac      = factor(N),
      rho_fac    = factor(target_rho),
      prior_fac  = factor(prior),
      method_fac = factor(method)
    ) %>%
    dplyr::filter(is.finite(log_metric))

  if (nrow(df) < 20) {
    warning("Too few observations for meta-regression: ", nrow(df))
    return(NULL)
  }

  # Fit OLS
  formula_str <- "log_metric ~ (N_fac + rho_fac + prior_fac + method_fac)^2"
  mod <- lm(as.formula(formula_str), data = df)

  # Cluster-robust SEs
  robust_vcov <- sandwich::vcovCL(
    mod,
    cluster = df$condition_id,
    type = "HC1"
  )
  robust_coefs <- lmtest::coeftest(mod, vcov. = robust_vcov)

  list(
    model        = mod,
    robust_coefs = robust_coefs,
    r_squared    = summary(mod)$r.squared,
    adj_r_sq     = summary(mod)$adj.r.squared,
    n_obs        = nrow(df),
    formula_used = formula_str
  )
}


# =============================================================================
# Plotting Functions
# =============================================================================

#' Plot Pilot Results: RMSE by Condition
#'
#' @param summary_df Summary data frame from pilot test.
#' @param save_path File path for saving figure.
#' @return ggplot object.
#' @export
plot_pilot_rmse <- function(summary_df, save_path = NULL) {

  df <- summary_df %>%
    dplyr::filter(!is.na(mean_rmse)) %>%
    dplyr::mutate(
      cond_label = sprintf("%s\n%s\nN=%d, rho=%.1f",
                           irt_model, latent_shape, N, target_rho),
      prior_label = dplyr::case_when(
        prior == "gaussian"   ~ "Gaussian",
        prior == "dp_broad"   ~ "DP-Broad",
        prior == "dp_focused" ~ "DP-Focused",
        TRUE ~ prior
      ),
      method_label = toupper(method)
    )

  p <- ggplot(df, aes(x = method_label, y = mean_rmse,
                       fill = prior_label)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    facet_wrap(~ cond_label, scales = "free_y", nrow = 1) +
    scale_fill_manual(
      values = c("Gaussian" = "#2166AC", "DP-Broad" = "#B2182B", "DP-Focused" = "#4DAF4A")
    ) +
    labs(
      title = "Pilot Results: Mean RMSE by Prior and Summary Method",
      x = "Summary Method",
      y = "Mean RMSE",
      fill = "Prior"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold"),
      strip.text = element_text(size = 8)
    )

  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 14, height = 5, dpi = 150)
  }

  p
}


#' Plot Pilot Results: KS Distance by Condition
#'
#' @param summary_df Summary data frame from pilot test.
#' @param save_path File path for saving figure.
#' @return ggplot object.
#' @export
plot_pilot_ks <- function(summary_df, save_path = NULL) {

  df <- summary_df %>%
    dplyr::filter(!is.na(mean_ks)) %>%
    dplyr::mutate(
      cond_label = sprintf("%s\n%s\nN=%d, rho=%.1f",
                           irt_model, latent_shape, N, target_rho),
      prior_label = dplyr::case_when(
        prior == "gaussian"   ~ "Gaussian",
        prior == "dp_broad"   ~ "DP-Broad",
        prior == "dp_focused" ~ "DP-Focused",
        TRUE ~ prior
      ),
      method_label = toupper(method)
    )

  p <- ggplot(df, aes(x = method_label, y = mean_ks,
                       fill = prior_label)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    facet_wrap(~ cond_label, scales = "free_y", nrow = 1) +
    scale_fill_manual(
      values = c("Gaussian" = "#2166AC", "DP-Broad" = "#B2182B", "DP-Focused" = "#4DAF4A")
    ) +
    labs(
      title = "Pilot Results: Mean KS Distance by Prior and Summary Method",
      x = "Summary Method",
      y = "Mean KS Distance",
      fill = "Prior"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold"),
      strip.text = element_text(size = 8)
    )

  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 14, height = 5, dpi = 150)
  }

  p
}
