cat("irw installed:", system.file(package = "irw") != "", "\n")
if (requireNamespace("irw", quietly = TRUE)) {
  cat("irw version:", as.character(packageVersion("irw")), "\n")
} else {
  cat("irw NOT available\n")
}

library(IRTsimrel)
cat("\neqc_calibrate signature:\n")
print(args(eqc_calibrate))

# Quick test: can we run with item_source = "irw"?
cat("\n--- Testing eqc_calibrate with item_source = 'irw' ---\n")
tryCatch({
  res <- eqc_calibrate(
    target_rho = 0.70,
    n_items = 25L,
    model = "rasch",
    latent_shape = "normal",
    item_source = "irw",
    reliability_metric = "info",
    M = 5000L,
    seed = 42L,
    verbose = TRUE
  )
  cat("SUCCESS! c* =", res$c_star, "\n")
  cat("achieved_rho =", res$achieved_rho, "\n")
  cat("beta range:", range(res$beta_vec), "\n")
  cat("lambda range:", range(res$lambda_scaled), "\n")
}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
})
