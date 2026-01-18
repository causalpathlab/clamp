#!/usr/bin/env Rscript
# Check if all required R packages and files are available

cat("Checking dependencies for clamp_binary()...\n\n")

# ============================================================================
# Check Required R Packages
# ============================================================================
cat("=== Checking R Packages ===\n")

required_packages <- c(
  "stats",       # Base R
  "Matrix",      # Matrix operations
  "matrixStats", # Column statistics
  "sandwich"     # Robust variance estimation (vcovHC)
)

missing_packages <- c()

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("✗ %s: NOT INSTALLED\n", pkg))
    missing_packages <- c(missing_packages, pkg)
  } else {
    cat(sprintf("✓ %s: OK\n", pkg))
  }
}

if (length(missing_packages) > 0) {
  cat("\nTo install missing packages, run:\n")
  cat(sprintf("  install.packages(c(%s))\n",
              paste(sprintf('"%s"', missing_packages), collapse = ", ")))
} else {
  cat("\nAll required packages are installed!\n")
}

# ============================================================================
# Check Required R Files
# ============================================================================
cat("\n=== Checking R Source Files ===\n")

required_files <- c(
  "../R/clamp_binary.R",
  "../R/clamp_update_each_effect_binary.R",
  "../R/ipw_single_effect_regression_binary.R",
  "../R/bootstrap_ipw_variance_binary.R",
  "../R/estimate_average_treatment_effect_binary.R",
  "../R/compute_colstats.R",
  "../R/initialize.R",
  "../R/elbo.R",
  "../R/robust_importance_weights.R",
  "../R/remove_abnormal_subjects.R",
  "../R/check_abnormal_subjects.R",
  "../R/sparse_multiplication.R",
  "../R/estimate_residual_variance.R",
  "../R/optimize_prior_variance.R",
  "../R/model_weighted_linear.R",
  "../R/clamp_utils.R",
  "../R/summary.clamp.R"
)

missing_files <- c()

for (file in required_files) {
  if (file.exists(file)) {
    cat(sprintf("✓ %s\n", file))
  } else {
    cat(sprintf("✗ %s: NOT FOUND\n", file))
    missing_files <- c(missing_files, file)
  }
}

if (length(missing_files) > 0) {
  cat("\nERROR: Some required files are missing!\n")
  cat("Make sure you're running this from the clamp/examples directory.\n")
} else {
  cat("\nAll required source files found!\n")
}

# ============================================================================
# Summary
# ============================================================================
cat("\n=== Summary ===\n")

if (length(missing_packages) == 0 && length(missing_files) == 0) {
  cat("✓ All dependencies satisfied!\n")
  cat("✓ You can now run: source('test_clamp_binary_simple.R')\n")
} else {
  cat("✗ Some dependencies are missing. Please resolve issues above.\n")
}
