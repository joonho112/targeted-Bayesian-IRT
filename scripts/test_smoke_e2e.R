#!/usr/bin/env Rscript
# =============================================================================
# test_smoke_e2e.R  —  End-to-End Smoke Test
# =============================================================================
# Runs the ENTIRE pipeline for 1 condition x 1 replication x 3 priors
# using SHORT MCMC chains (500 iter) to verify everything connects.
#
# Expected runtime: ~3-5 minutes
# =============================================================================

cat("\n====================================================\n")
cat("  END-TO-END SMOKE TEST\n")
cat("  1 condition, 1 rep, 3 priors (short MCMC)\n")
cat("====================================================\n\n")

library(nimble)       # MUST load before DPMirt
library(DPMirt)
library(DPprior)
library(IRTsimrel)
library(tidyverse)
library(yaml)
library(digest)

source("R/00_config.R")
source("R/01_calibrate.R")
source("R/02_generate.R")
source("R/03_fit.R")

cfg <- load_config()

# Override MCMC settings for speed
cfg$mcmc$niter   <- 500L
cfg$mcmc$nburnin <- 100L
cfg$posterior_storage$save_posteriors <- FALSE

# =============================================================================
# Step 1: Pick one condition
# =============================================================================
design <- create_design_matrix(cfg)
# Choose: Rasch, normal, N=100, rho=0.70 (condition ~middle of the grid)
test_row <- design %>%
  filter(irt_model == "rasch", latent_shape == "normal",
         N == 100, target_rho == 0.70) %>%
  slice(1)
cat(sprintf("Test condition: %s, %s, N=%d, rho=%.2f (id=%d)\n\n",
    test_row$irt_model, test_row$latent_shape, test_row$N,
    test_row$target_rho, test_row$condition_id))

# =============================================================================
# Step 2: Calibrate
# =============================================================================
cat("--- Step 2: EQC Calibration ---\n")
t0 <- Sys.time()
calib_results <- run_all_calibrations(cfg)
eqc_result <- lookup_calibration(test_row, calib_results)
cat(sprintf("  c* = %.4f, achieved rho = %.4f\n",
    eqc_result$c_star, eqc_result$achieved_rho))
cat(sprintf("  Time: %.1f sec\n\n", difftime(Sys.time(), t0, units="secs")))

# =============================================================================
# Step 3: Generate one dataset
# =============================================================================
cat("--- Step 3: Data Generation ---\n")
rep_seed <- derive_seed(cfg$seeds$master_seed, test_row$condition_id, rep_id = 1L)
dataset <- generate_one_dataset(
  eqc_result   = eqc_result,
  N            = test_row$N,
  latent_shape = test_row$latent_shape,
  rep_seed     = rep_seed
)
cat(sprintf("  Response matrix: %d x %d\n", dataset$N, dataset$I))
cat(sprintf("  theta range: [%.3f, %.3f], mean=%.3f, sd=%.3f\n",
    min(dataset$theta), max(dataset$theta),
    mean(dataset$theta), sd(dataset$theta)))
cat(sprintf("  Item p-values range: [%.2f, %.2f]\n",
    min(colMeans(dataset$response_matrix)),
    max(colMeans(dataset$response_matrix))))

# =============================================================================
# Step 4: Build alpha prior cache
# =============================================================================
cat("\n--- Step 4: Alpha Prior Cache ---\n")
alpha_cache <- build_alpha_cache(cfg)
cat(sprintf("  Cached %d DP prior specifications\n", length(alpha_cache)))
for (key in names(alpha_cache)) {
  ab <- alpha_cache[[key]]
  cat(sprintf("    %s: Gamma(a=%.4f, b=%.4f), E[alpha]=%.4f\n",
      key, ab["a"], ab["b"], ab["a"]/ab["b"]))
}

# =============================================================================
# Step 5: Fit all 3 priors
# =============================================================================
cat("\n--- Step 5: Model Fitting (3 priors) ---\n")
prior_specs <- get_prior_specs(cfg, test_row$N)

all_results <- list()
for (prior_name in names(prior_specs)) {
  pspec <- prior_specs[[prior_name]]
  cat(sprintf("\n  Fitting %s prior...\n", pspec$label))
  t0 <- Sys.time()

  model_seed <- derive_seed(cfg$seeds$master_seed, test_row$condition_id,
                            rep_id = 1L, prior_index = pspec$prior_index)

  alpha_prior <- NULL
  if (pspec$type == "dpm") {
    cache_key <- sprintf("%s_N%d", prior_name, test_row$N)
    alpha_prior <- alpha_cache[[cache_key]]
  }

  result <- fit_one_model(
    dataset     = dataset,
    irt_model   = test_row$irt_model,
    prior_spec  = pspec,
    alpha_prior = alpha_prior,
    model_seed  = model_seed,
    config      = cfg
  )

  elapsed <- difftime(Sys.time(), t0, units = "secs")
  cat(sprintf("    Status: %s (%.1f sec)\n", result$fit_summary$status, elapsed))

  if (result$fit_summary$status == "success") {
    # Show theta estimates
    theta_pm <- result$estimates$theta$theta_pm
    cat(sprintf("    theta_pm range: [%.3f, %.3f]\n", min(theta_pm), max(theta_pm)))

    # Show losses
    theta_loss <- result$losses[result$losses$parameter == "theta", ]
    for (j in 1:nrow(theta_loss)) {
      cat(sprintf("    %s: RMSE=%.4f, MSELR=%.6f, KS=%.4f\n",
          toupper(theta_loss$method[j]),
          sqrt(theta_loss$msel[j]),
          theta_loss$mselr[j],
          theta_loss$ks[j]))
    }
  }

  all_results[[prior_name]] <- result
}

# =============================================================================
# Step 6: Summary & Sanity Checks
# =============================================================================
cat("\n====================================================\n")
cat("  SANITY CHECKS\n")
cat("====================================================\n\n")

# Check 1: All fits succeeded
all_ok <- all(sapply(all_results, function(r) r$fit_summary$status == "success"))
cat(sprintf("  [%s] All 3 fits succeeded\n", ifelse(all_ok, "PASS", "FAIL")))

# Check 2: PM should have lowest RMSE (it's the Bayes optimal for MSE)
if (all_ok) {
  for (pname in names(all_results)) {
    tl <- all_results[[pname]]$losses
    tl <- tl[tl$parameter == "theta", ]
    rmse_pm <- sqrt(tl$msel[tl$method == "pm"])
    rmse_cb <- sqrt(tl$msel[tl$method == "cb"])
    rmse_gr <- sqrt(tl$msel[tl$method == "gr"])
    pm_best <- rmse_pm <= rmse_cb && rmse_pm <= rmse_gr
    cat(sprintf("  [%s] PM has lowest RMSE for %s (PM=%.4f, CB=%.4f, GR=%.4f)\n",
        ifelse(pm_best, "PASS", "NOTE"), pname, rmse_pm, rmse_cb, rmse_gr))
  }

  # Check 3: GR should have best KS (it's designed for distributional recovery)
  for (pname in names(all_results)) {
    tl <- all_results[[pname]]$losses
    tl <- tl[tl$parameter == "theta", ]
    ks_pm <- tl$ks[tl$method == "pm"]
    ks_gr <- tl$ks[tl$method == "gr"]
    gr_better_ks <- ks_gr <= ks_pm
    cat(sprintf("  [%s] GR has better KS than PM for %s (GR=%.4f, PM=%.4f)\n",
        ifelse(gr_better_ks, "PASS", "NOTE"), pname, ks_gr, ks_pm))
  }
}

cat("\n====================================================\n")
cat("  END-TO-END SMOKE TEST COMPLETE\n")
cat("====================================================\n")
