# =============================================================================
# 04_generate_figures.R
# =============================================================================
# Phase 4: Generate publication-quality figures.
#
# Figures produced:
#   1. Loss metric boxplots across priors and summary methods
#   2. Meta-regression coefficient plots
#   3. Reliability vs. loss line plots
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


# =============================================================================
# 1. Load Post-Processed Data
# =============================================================================

results_summary <- readRDS(
  file.path(config$paths$results_dir, "results_summary_mean.rds")
)

meta_results <- readRDS(
  file.path(config$paths$results_dir, "meta_regression_results.rds")
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
# Figure 1: Loss Metric Boxplots
# =============================================================================
#
# Boxplots of loss metrics across priors and summary methods,
# faceted by latent_shape and target_rho.

plot_loss_boxplots <- function(results_df,
                               metric_name = "rmse",
                               model_fix   = "rasch") {

  df_sub <- results_df |>
    dplyr::filter(
      irt_model == model_fix,
      status    == "success"
    )

  if (!(metric_name %in% names(df_sub))) {
    warning("Metric '", metric_name, "' not found in results.")
    return(NULL)
  }

  p <- ggplot(df_sub, aes(
    x    = method,
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
      subtitle = paste0("Model: ", toupper(model_fix))
    ) +
    theme_pub

  p
}


# Generate boxplot figures
for (metric in config$metrics) {
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
# Figure 2: Meta-Regression Coefficient Summary
# =============================================================================
#
# For each metric's meta-regression, extract the top coefficients and
# plot their magnitudes.

plot_meta_coefs <- function(meta_result, metric_name, n_top = 15) {

  if (is.null(meta_result)) return(NULL)

  coef_tbl <- broom::tidy(meta_result$robust_coefs)
  coef_tbl <- coef_tbl |>
    dplyr::filter(term != "(Intercept)") |>
    dplyr::arrange(dplyr::desc(abs(estimate))) |>
    dplyr::slice_head(n = n_top) |>
    dplyr::mutate(
      term = factor(term, levels = rev(term)),
      direction = ifelse(estimate > 0, "Increase", "Decrease")
    )

  p <- ggplot(coef_tbl, aes(x = estimate, y = term, fill = direction)) +
    geom_col() +
    geom_errorbarh(
      aes(xmin = estimate - 1.96 * std.error,
          xmax = estimate + 1.96 * std.error),
      height = 0.3
    ) +
    scale_fill_manual(
      values = c("Increase" = "firebrick", "Decrease" = "steelblue")
    ) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(
      x     = "Coefficient (log scale)",
      y     = NULL,
      fill  = "Direction",
      title = paste0("Meta-Regression: Top ", n_top, " Coefficients"),
      subtitle = paste0("Metric: ", toupper(metric_name),
                        " | R^2 = ", sprintf("%.3f", meta_result$r_squared))
    ) +
    theme_pub +
    theme(legend.position = "right")

  p
}


for (metric in names(meta_results)) {
  p <- plot_meta_coefs(meta_results[[metric]], metric)
  if (!is.null(p)) {
    fname <- sprintf("meta_coefs_%s.pdf", metric)
    ggsave(
      file.path(config$paths$figures_dir, fname),
      p, width = 10, height = 7
    )
  }
}

log_msg("Meta-regression coefficient figures saved.", config = config)


# =============================================================================
# Figure 3: Loss vs. Reliability (rho)
# =============================================================================
#
# Line plots showing how each loss metric changes as a function of
# target reliability, by prior and summary method.

plot_loss_vs_rho <- function(summary_df,
                             metric_name = "rmse",
                             model_fix   = "rasch") {

  # Check for the expected column naming pattern
  # summarize_results() produces columns like rmse_mean, rmse_sd, etc.
  metric_mean_col <- paste0(metric_name, "_mean")

  if (!(metric_mean_col %in% names(summary_df))) {
    # Try bare metric name in case column naming differs
    if (metric_name %in% names(summary_df)) {
      metric_mean_col <- metric_name
    } else {
      warning("Column '", metric_mean_col, "' not found. Skipping.")
      return(NULL)
    }
  }

  df_sub <- summary_df |>
    dplyr::filter(irt_model == model_fix)

  if (nrow(df_sub) == 0) return(NULL)

  p <- ggplot(df_sub, aes(
    x        = target_rho,
    y        = .data[[metric_mean_col]],
    colour   = prior,
    linetype = method
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
      subtitle = paste0("Model: ", toupper(model_fix))
    ) +
    theme_pub

  p
}


for (metric in config$metrics) {
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
