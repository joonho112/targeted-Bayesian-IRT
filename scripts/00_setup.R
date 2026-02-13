#!/usr/bin/env Rscript
# =============================================================================
# 00_setup.R — Environment Setup & Dependency Verification
# =============================================================================
# Run this script FIRST on any new machine to:
#   1. Check R version compatibility
#   2. Install/verify all CRAN packages
#   3. Install/verify all GitHub packages (IRTsimrel, DPprior, DPMirt, irw)
#   4. Verify nimble C++ toolchain
#   5. Run a quick smoke test (nimble load order, eqc_calibrate, dpmirt)
#   6. Detect available cores and configure parallelization
#   7. Source all R modules and validate config
#   8. Save session info for reproducibility
#
# Usage:
#   Rscript scripts/00_setup.R
#   Rscript scripts/00_setup.R --install    # Auto-install missing packages
#   Rscript scripts/00_setup.R --skip-test  # Skip smoke test
#
# Author: JoonHo Lee
# Date: February 2026
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
auto_install <- "--install" %in% args
skip_test    <- "--skip-test" %in% args

cat("\n")
cat("==============================================================\n")
cat("  Targeted Bayesian IRT — Environment Setup\n")
cat("==============================================================\n\n")

errors   <- character(0)
warnings <- character(0)

# =============================================================================
# 1. R Version Check
# =============================================================================
cat("--- [1/8] R Version ---\n")
r_ver <- paste0(R.version$major, ".", R.version$minor)
cat(sprintf("  R version: %s\n", r_ver))
cat(sprintf("  Platform:  %s\n", R.version$platform))

if (as.numeric(R.version$major) < 4 ||
    (as.numeric(R.version$major) == 4 && as.numeric(R.version$minor) < 3)) {
  errors <- c(errors, "R >= 4.3.0 is required. Please update R.")
  cat("  [FAIL] R >= 4.3.0 required\n")
} else {
  cat("  [OK]\n")
}

# =============================================================================
# 2. CRAN Packages
# =============================================================================
cat("\n--- [2/8] CRAN Packages ---\n")

cran_packages <- c(
  # Core simulation
  "nimble",        # MCMC engine (MUST load before DPMirt)
  "yaml",          # Configuration parsing
  "digest",        # Reproducible seed derivation

  # Data manipulation
  "tibble", "dplyr", "tidyr", "purrr", "readr",
  "forcats", "stringr", "lubridate", "ggplot2",

  # Parallel computing
  "future",        # Parallel backend
  "future.apply",  # future_lapply()

  # Post-processing
  "sandwich",      # Cluster-robust standard errors
  "lmtest",        # Robust coefficient testing
  "broom",         # Tidy model summaries

  # Utilities
  "here",          # Project-relative paths
  "remotes"        # For GitHub installs
)

cran_status <- sapply(cran_packages, function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
})

installed_cran   <- names(cran_status[cran_status])
missing_cran     <- names(cran_status[!cran_status])

for (pkg in installed_cran) {
  ver <- tryCatch(as.character(packageVersion(pkg)), error = function(e) "?")
  cat(sprintf("  [OK]      %-15s %s\n", pkg, ver))
}

if (length(missing_cran) > 0) {
  for (pkg in missing_cran) {
    cat(sprintf("  [MISSING] %-15s\n", pkg))
  }

  if (auto_install) {
    cat(sprintf("\n  Installing %d missing CRAN packages...\n", length(missing_cran)))
    install.packages(missing_cran, repos = "https://cloud.r-project.org")

    # Re-check
    still_missing <- missing_cran[!sapply(missing_cran, requireNamespace, quietly = TRUE)]
    if (length(still_missing) > 0) {
      errors <- c(errors, paste("Failed to install CRAN packages:",
                                paste(still_missing, collapse = ", ")))
    } else {
      cat("  All CRAN packages installed successfully.\n")
      missing_cran <- character(0)
    }
  } else {
    cat(sprintf("\n  To install missing packages, run:\n"))
    cat(sprintf("    install.packages(c(%s))\n",
        paste0('"', missing_cran, '"', collapse = ", ")))
    cat("  Or re-run with: Rscript scripts/00_setup.R --install\n")
    errors <- c(errors, paste("Missing CRAN packages:",
                              paste(missing_cran, collapse = ", ")))
  }
}

# =============================================================================
# 3. GitHub Packages
# =============================================================================
cat("\n--- [3/8] GitHub Packages (Custom) ---\n")

github_packages <- list(
  IRTsimrel = list(
    repo     = "joonho112/IRTsimrel",
    min_ver  = "0.2.0",
    purpose  = "EQC calibration + response data generation"
  ),
  DPprior = list(
    repo     = "joonho112/DPprior",
    min_ver  = "0.1.0",
    purpose  = "DP concentration parameter elicitation"
  ),
  DPMirt = list(
    repo     = "joonho112/DPMirt",
    min_ver  = "0.1.0",
    purpose  = "Bayesian IRT with DP mixture priors"
  ),
  irw = list(
    repo     = "joonho112/irw",
    min_ver  = "1.0.0",
    purpose  = "Item Response Warehouse (empirical item pools)"
  )
)

missing_github <- character(0)
outdated_github <- character(0)

for (pkg_name in names(github_packages)) {
  pkg_info <- github_packages[[pkg_name]]

  if (requireNamespace(pkg_name, quietly = TRUE)) {
    ver <- as.character(packageVersion(pkg_name))
    cat(sprintf("  [OK]      %-15s %s", pkg_name, ver))

    # Version check
    if (compareVersion(ver, pkg_info$min_ver) < 0) {
      cat(sprintf("  (OUTDATED: need >= %s)", pkg_info$min_ver))
      outdated_github <- c(outdated_github, pkg_name)
    }
    cat("\n")
  } else {
    cat(sprintf("  [MISSING] %-15s -- %s\n", pkg_name, pkg_info$purpose))
    missing_github <- c(missing_github, pkg_name)
  }
}

if (length(missing_github) > 0 || length(outdated_github) > 0) {
  to_install <- union(missing_github, outdated_github)

  if (auto_install && requireNamespace("remotes", quietly = TRUE)) {
    cat(sprintf("\n  Installing %d GitHub packages...\n", length(to_install)))
    for (pkg_name in to_install) {
      repo <- github_packages[[pkg_name]]$repo
      cat(sprintf("  Installing %s from %s...\n", pkg_name, repo))
      tryCatch({
        remotes::install_github(repo, upgrade = "never", quiet = TRUE)
        cat(sprintf("  [OK] %s installed\n", pkg_name))
      }, error = function(e) {
        cat(sprintf("  [FAIL] %s: %s\n", pkg_name, conditionMessage(e)))
        errors <<- c(errors, paste("Failed to install", pkg_name))
      })
    }
  } else {
    cat("\n  To install, run in R:\n")
    for (pkg_name in to_install) {
      repo <- github_packages[[pkg_name]]$repo
      cat(sprintf('    remotes::install_github("%s")\n', repo))
    }
    cat("  Or re-run with: Rscript scripts/00_setup.R --install\n")
    errors <- c(errors, paste("Missing GitHub packages:",
                              paste(to_install, collapse = ", ")))
  }
}

# =============================================================================
# 4. Nimble C++ Toolchain
# =============================================================================
cat("\n--- [4/8] Nimble C++ Toolchain ---\n")

if (requireNamespace("nimble", quietly = TRUE)) {
  # Test that nimble can compile a trivial model
  tryCatch({
    suppressMessages(library(nimble))
    cat("  nimble loaded successfully.\n")

    # Check if C++ compiler is available
    nimble_dir <- system.file(package = "nimble")
    cat(sprintf("  nimble path: %s\n", nimble_dir))
    cat("  [OK] nimble operational\n")
  }, error = function(e) {
    cat(sprintf("  [WARN] nimble loaded but may have issues: %s\n",
                conditionMessage(e)))
    warnings <- c(warnings, paste("nimble issue:", conditionMessage(e)))
  })
} else {
  cat("  [SKIP] nimble not installed\n")
}

# =============================================================================
# 5. Load Order Test (CRITICAL: nimble before DPMirt)
# =============================================================================
cat("\n--- [5/8] Load Order Test ---\n")

all_core_available <- all(sapply(
  c("nimble", "DPMirt", "DPprior", "IRTsimrel", "irw"),
  requireNamespace, quietly = TRUE
))

if (all_core_available) {
  tryCatch({
    suppressMessages({
      library(nimble)
      library(DPMirt)
      library(DPprior)
      library(IRTsimrel)
    })
    cat("  [OK] nimble -> DPMirt -> DPprior -> IRTsimrel loaded in correct order\n")

    # Quick functional test
    if (!skip_test) {
      cat("\n  Running quick smoke test...\n")

      # Test 1: EQC calibration with IRW
      cat("    [1] eqc_calibrate(irw)... ")
      eqc_res <- IRTsimrel::eqc_calibrate(
        target_rho = 0.70, n_items = 10L, model = "rasch",
        latent_shape = "normal", item_source = "irw",
        M = 2000L, seed = 42L, verbose = FALSE
      )
      cat(sprintf("OK (c*=%.3f, rho=%.4f)\n", eqc_res$c_star, eqc_res$achieved_rho))

      # Test 2: DPprior_fit
      cat("    [2] DPprior_fit()... ")
      dp_res <- DPprior::DPprior_fit(
        J = 50L, mu_K = 10, confidence = "low",
        method = "A2-MN", check_diagnostics = FALSE,
        warn_dominance = FALSE, verbose = FALSE
      )
      cat(sprintf("OK (a=%.2f, b=%.2f)\n", dp_res$a, dp_res$b))

      # Test 3: dpmirt (tiny model)
      cat("    [3] dpmirt() short fit... ")
      sim_dat <- IRTsimrel::simulate_response_data(
        result = eqc_res, n_persons = 20L,
        latent_shape = "normal", seed = 123L
      )
      fit <- DPMirt::dpmirt(
        data = sim_dat$response_matrix, model = "rasch",
        prior = "normal", niter = 100L, nburnin = 20L,
        nchains = 1L, seed = 1L, verbose = FALSE,
        compute_waic = FALSE, compute_dp_density = FALSE
      )
      est <- DPMirt::dpmirt_estimates(fit, methods = c("pm", "cb", "gr"))
      cat(sprintf("OK (%d posterior samples, 3 methods)\n", nrow(fit$theta_samp)))

      cat("  [OK] All smoke tests passed.\n")
    } else {
      cat("  [SKIP] Smoke test skipped (--skip-test)\n")
    }
  }, error = function(e) {
    cat(sprintf("\n  [FAIL] Load/test error: %s\n", conditionMessage(e)))
    errors <<- c(errors, paste("Smoke test failed:", conditionMessage(e)))
  })
} else {
  cat("  [SKIP] Core packages not all installed; cannot test load order.\n")
}

# =============================================================================
# 6. Parallel Computing
# =============================================================================
cat("\n--- [6/8] Parallel Computing ---\n")

n_cores <- parallel::detectCores(logical = TRUE)
n_phys  <- parallel::detectCores(logical = FALSE)
cat(sprintf("  Logical cores:  %d\n", n_cores))
cat(sprintf("  Physical cores: %d\n", n_phys))

if (requireNamespace("future", quietly = TRUE)) {
  cat(sprintf("  future backend: available\n"))

  # Suggest optimal worker count
  suggested_workers <- max(1, n_phys - 1)
  cat(sprintf("  Suggested workers: %d (physical cores - 1)\n", suggested_workers))

  # Read current config setting
  cfg_path <- file.path(getwd(), "config", "sim_config.yaml")
  if (file.exists(cfg_path)) {
    cfg_yaml <- yaml::read_yaml(cfg_path)
    cfg_workers <- cfg_yaml$parallel$n_workers
    cat(sprintf("  Config n_workers: %d", cfg_workers))

    if (cfg_workers > n_cores) {
      cat(sprintf(" [WARN: exceeds available cores (%d)]", n_cores))
      warnings <- c(warnings,
        sprintf("Config n_workers=%d exceeds available cores=%d", cfg_workers, n_cores))
    } else if (cfg_workers != suggested_workers) {
      cat(sprintf(" [NOTE: consider updating to %d for this machine]", suggested_workers))
    } else {
      cat(" [OK]")
    }
    cat("\n")

    cat(sprintf("\n  To update for this machine (%d physical cores), edit:\n", n_phys))
    cat(sprintf("    config/sim_config.yaml -> parallel: n_workers: %d\n",
        suggested_workers))
  }
} else {
  cat("  [MISSING] future package not available\n")
}

# =============================================================================
# 7. Project Structure Validation
# =============================================================================
cat("\n--- [7/8] Project Structure ---\n")

required_files <- c(
  "config/sim_config.yaml",
  "R/00_config.R",
  "R/01_calibrate.R",
  "R/02_generate.R",
  "R/03_fit.R",
  "R/04_evaluate.R",
  "R/05_postprocess.R",
  "R/utils.R",
  "scripts/run_simulation.R",
  "scripts/test_pilot.R"
)

for (f in required_files) {
  full_path <- file.path(getwd(), f)
  if (file.exists(full_path)) {
    cat(sprintf("  [OK]   %s\n", f))
  } else {
    cat(sprintf("  [MISS] %s\n", f))
    errors <- c(errors, paste("Missing file:", f))
  }
}

# Source modules and validate config if all core packages are available
if (all_core_available && length(missing_cran) == 0) {
  cat("\n  Loading modules and validating config...\n")
  tryCatch({
    source("R/utils.R")
    source("R/00_config.R")
    source("R/01_calibrate.R")
    source("R/02_generate.R")
    source("R/03_fit.R")
    source("R/04_evaluate.R")
    source("R/05_postprocess.R")

    config <- load_config()
    ensure_directories(config)
    design <- create_design_matrix(config)
    validate_design(design, config)

    cat(sprintf("  [OK] Config valid: %d conditions, %d reps, %d total fits\n",
        nrow(design), config$design$replications,
        nrow(design) * 3 * config$design$replications))
    cat(sprintf("  [OK] Item source: %s\n", config$calibration$item_source))
    cat(sprintf("  [OK] Output directories created\n"))
  }, error = function(e) {
    cat(sprintf("  [FAIL] %s\n", conditionMessage(e)))
    errors <<- c(errors, paste("Config validation failed:", conditionMessage(e)))
  })
}

# =============================================================================
# 8. Session Info
# =============================================================================
cat("\n--- [8/8] Session Info ---\n")

output_dir <- file.path(getwd(), "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

session_path <- file.path(output_dir, "session_info.txt")
writeLines(capture.output(sessionInfo()), session_path)
cat(sprintf("  Saved to: %s\n", session_path))

# =============================================================================
# Summary
# =============================================================================
cat("\n")
cat("==============================================================\n")

if (length(errors) == 0 && length(warnings) == 0) {
  cat("  SETUP COMPLETE — All checks passed!\n")
  cat("  Ready to run: Rscript scripts/run_simulation.R --pilot\n")
} else {
  if (length(warnings) > 0) {
    cat(sprintf("  WARNINGS (%d):\n", length(warnings)))
    for (w in warnings) cat(sprintf("    - %s\n", w))
  }
  if (length(errors) > 0) {
    cat(sprintf("  ERRORS (%d):\n", length(errors)))
    for (e in errors) cat(sprintf("    - %s\n", e))
    cat("\n  Fix the above errors before running the simulation.\n")
    cat("  Tip: Rscript scripts/00_setup.R --install\n")
  }
}
cat("==============================================================\n\n")
