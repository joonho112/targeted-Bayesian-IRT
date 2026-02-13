# =============================================================================
# 04_generate_figures.R
# =============================================================================
# Phase 4: Generate publication-quality figures.
#
# Figures produced:
#   1. Scaled EDF panels (estimated vs. true density per condition)
#   2. Loss metric boxplots across priors and summary methods
#   3. Meta-regression heatmaps of marginal effects
#   4. Reliability vs. loss scatter plots
#
# Prerequisites:
#   - 03_run_postprocessing.R has been run (post-processed outputs on disk).
#
# Author: JoonHo Lee
# Date: February 2026
# =============================================================================


# --- Source setup if not already loaded ---
if (!exists("config")) {
  source(file.path(here::here(), "scripts", "00_setup.R"))
}

library(ggplot2)
library(patchwork)


# =============================================================================
# 1. Load Post-Processed Data
# =============================================================================

results_summary <- readRDS(
  file.path(config$paths$results_dir, "results_summary_mean.rds")
)

edf_results <- readRDS(
  file.path(config$paths$results_dir, "scaled_edf.rds")
)

meta_reg_results <- readRDS(
  file.path(config$paths$results_dir, "meta_regression_results.rds")
)

marginal_effects <- readRDS(
  file.path(config$paths$results_dir, "marginal_effects.rds")
)

# Also load full results for boxplots
results_all <- readRDS(
  file.path(config$paths$results_dir, "simulation_results_compiled.rds")
)


# =============================================================================
# Shared Theme
# =============================================================================

theme_pub <- theme_minimal(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom"
  )


# =============================================================================
# Figure 1: Scaled EDF Panels
# =============================================================================
#
# For a selected condition (e.g., bimodal, N=200, rho=0.8), show
# the scaled EDF across all (prior x summary_method) combinations.

plot_scaled_edf <- function(edf_df,
                            shape_fix   = "bimodal",
                            N_fix       = 200,
                            rho_fix     = 0.8,
                            model_fix   = "rasch",
                            abs_limit   = 5) {

  df_sub <- edf_df |>
    dplyr::filter(
      latent_shape == shape_fix,
      N            == N_fix,
      target_rho   == rho_fix,
      irt_model    == model_fix,
      abs(bin_mid) <= abs_limit
    )

  if (nrow(df_sub) == 0) {
    warning("No data for the specified condition.")
    return(NULL)
  }

  subtitle_lab <- sprintf(
    "G = %s, N = %d, rho = %.1f, Model = %s",
    shape_fix, N_fix, rho_fix, toupper(model_fix)
  )

  p <- ggplot(df_sub) +
    geom_rect(
      aes(xmin = bin_start, xmax = bin_end,
          ymin = 0, ymax = density_est),
      fill = "limegreen", colour = "grey80", linewidth = 0.1
    ) +
    geom_line(
      aes(x = bin_mid, y = density_true_sp),
      linewidth = 0.6, colour = "black"
    ) +
    geom_line(
      aes(x = bin_mid, y = density_true_fs),
      linewidth = 0.6, linetype = "dashed", colour = "black"
    ) +
    facet_grid(
      rows = vars(prior),
      cols = vars(summary_method),
      labeller = labeller(.default = label_value)
    ) +
    labs(
      x        = "Estimated Latent Parameter",
      y        = "Density",
      title    = "Scaled Empirical Distribution Function vs. True Distribution",
      subtitle = subtitle_lab,
      caption  = "Solid line: super-population; Dashed: finite-sample true"
    ) +
    theme_pub

  p
}


# Generate EDF figures for key conditions
edf_conditions <- tidyr::expand_grid(
  latent_shape = config$design$latent_shapes,
  N            = c(100, 500),
  target_rho   = c(0.7, 0.9),
  irt_model    = "rasch"
)

for (i in seq_len(nrow(edf_conditions))) {
  ec <- edf_conditions[i, ]

  p <- plot_scaled_edf(
    edf_results,
    shape_fix = ec$latent_shape,
    N_fix     = ec$N,
    rho_fix   = ec$target_rho,
    model_fix = ec$irt_model
  )

  if (!is.null(p)) {
    fname <- sprintf(
      "edf_%s_N%d_rho%.1f_%s.pdf",
      ec$latent_shape, ec$N, ec$target_rho, ec$irt_model
    )
    ggsave(
      file.path(config$paths$figures_dir, fname),
      p, width = 10, height = 8
    )
  }
}

log_msg("Scaled EDF figures saved.", config = config)


# =============================================================================
# Figure 2: Loss Metric Boxplots
# =============================================================================
#
# Boxplots of loss metrics across priors and summary methods,
# faceted by latent_shape and target_rho.

plot_loss_boxplots <- function(results_df,
                               metric_name = "rmse",
                               param_fix   = "theta",
                               model_fix   = "rasch") {

  df_sub <- results_df |>
    dplyr::filter(
      parameter == param_fix,
      irt_model == model_fix,
      status    == "success"
    )

  if (!(metric_name %in% names(df_sub))) {
    warning("Metric '", metric_name, "' not found in results.")
    return(NULL)
  }

  p <- ggplot(df_sub, aes(
    x    = summary_method,
    y    = .data[[metric_name]],
    fill = prior
  )) +
    geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
    facet_grid(
      rows = vars(latent_shape),
      cols = vars(target_rho),
      scales = "free_y",
      labeller = label_both
    ) +
    scale_y_log10() +
    labs(
      x     = "Summary Method",
      y     = paste0(toupper(metric_name), " (log scale)"),
      fill  = "Prior",
      title = paste0(toupper(metric_name), " by Prior and Summary Method"),
      subtitle = paste0("Parameter: ", param_fix,
                        ", Model: ", toupper(model_fix))
    ) +
    theme_pub

  p
}


# Generate boxplot figures
for (metric in c("rmse", "mselp", "ks")) {
  for (model in config$design$irt_models) {
    p <- plot_loss_boxplots(results_all, metric_name = metric,
                            model_fix = model)
    if (!is.null(p)) {
      fname <- sprintf("boxplot_%s_%s.pdf", metric, model)
      ggsave(
        file.path(config$paths$figures_dir, fname),
        p, width = 14, height = 10
      )
    }
  }
}

log_msg("Loss boxplot figures saved.", config = config)


# =============================================================================
# Figure 3: Meta-Regression Heatmap of Single Factor Effects
# =============================================================================
#
# Replicates the old codebase's plot_single_effects():
# Heatmap of multiplicative change in each metric when a single factor
# changes from its reference level.

plot_marginal_heatmap <- function(marginal_df,
                                  shape_fix = "bimodal",
                                  model_fix = "rasch",
                                  multiplicative = TRUE) {

  df_sub <- marginal_df |>
    dplyr::filter(
      latent_shape == shape_fix,
      irt_model    == model_fix
    )

  if (nrow(df_sub) == 0) return(NULL)

  df_sub$cell <- if (multiplicative) df_sub$exp_estimate else df_sub$estimate
  midpoint    <- if (multiplicative) 1 else 0
  fill_label  <- if (multiplicative) "Multiplicative\nchange" else "Log change"

  p <- ggplot(df_sub, aes(x = metric, y = term, fill = cell)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", cell)), size = 3) +
    scale_fill_gradient2(
      low = "steelblue", mid = "white", high = "firebrick",
      midpoint = midpoint
    ) +
    labs(
      x       = "Metric",
      y       = "Factor",
      fill    = fill_label,
      title   = "Meta-Regression: Single Factor Effects",
      caption = paste0("G = ", shape_fix, ", Model = ", toupper(model_fix))
    ) +
    theme_pub +
    theme(legend.position = "right")

  p
}


for (shape in config$design$latent_shapes) {
  for (model in config$design$irt_models) {
    p <- plot_marginal_heatmap(marginal_effects,
                               shape_fix = shape,
                               model_fix = model)
    if (!is.null(p)) {
      fname <- sprintf("heatmap_marginal_%s_%s.pdf", shape, model)
      ggsave(
        file.path(config$paths$figures_dir, fname),
        p, width = 8, height = 6
      )
    }
  }
}

log_msg("Meta-regression heatmap figures saved.", config = config)


# =============================================================================
# Figure 4: Loss vs. Reliability (rho)
# =============================================================================
#
# Line plots showing how each loss metric changes as a function of
# target reliability, by prior and summary method.

plot_loss_vs_rho <- function(summary_df,
                             metric_name = "rmse",
                             param_fix   = "theta",
                             model_fix   = "rasch") {

  metric_mean_col <- paste0(metric_name, "_mean")

  if (!(metric_mean_col %in% names(summary_df))) return(NULL)

  df_sub <- summary_df |>
    dplyr::filter(
      parameter == param_fix,
      irt_model == model_fix
    )

  p <- ggplot(df_sub, aes(
    x      = target_rho,
    y      = .data[[metric_mean_col]],
    colour = prior,
    linetype = summary_method
  )) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    facet_grid(
      rows   = vars(latent_shape),
      cols   = vars(N),
      scales = "free_y",
      labeller = label_both
    ) +
    scale_y_log10() +
    labs(
      x        = "Target Reliability (rho)",
      y        = paste0("Mean ", toupper(metric_name), " (log scale)"),
      colour   = "Prior",
      linetype = "Summary",
      title    = paste0(toupper(metric_name), " vs. Reliability"),
      subtitle = paste0("Parameter: ", param_fix,
                        ", Model: ", toupper(model_fix))
    ) +
    theme_pub

  p
}


for (metric in c("rmse", "mselp", "ks")) {
  for (model in config$design$irt_models) {
    p <- plot_loss_vs_rho(results_summary, metric_name = metric,
                          model_fix = model)
    if (!is.null(p)) {
      fname <- sprintf("loss_vs_rho_%s_%s.pdf", metric, model)
      ggsave(
        file.path(config$paths$figures_dir, fname),
        p, width = 14, height = 10
      )
    }
  }
}

log_msg("Loss vs. rho figures saved.", config = config)


# =============================================================================
# Summary
# =============================================================================

figure_files <- list.files(config$paths$figures_dir, pattern = "\\.pdf$")
cat(sprintf("\n=== Figure Generation Complete ===\n"))
cat(sprintf("  %d figures saved to: %s\n",
            length(figure_files), config$paths$figures_dir))
