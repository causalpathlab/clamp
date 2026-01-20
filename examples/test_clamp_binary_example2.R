#!/usr/bin/env Rscript
# Test example 2 for clamp_binary() function
# Tests with larger dataset and more causal variables
# High-dimensional cases
# Not converging.

# Load required packages
library(Matrix)
library(matrixStats)
library(sandwich)

set.seed(789)

# ============================================================================
# 1. Simulate Data, a high-dimensional case
# ============================================================================

n <- 1500  # number of samples
p <- 50    # number of treatments

# Generate 1-dimensional confounder
Z <- rnorm(n, mean = 0, sd = 1)

# ============================================================================
# 2. Generate Binary Treatments (confounded by Z)
# ============================================================================
# Each treatment X[,j] is generated from logistic regression with confounder Z
# P(X[,j] = 1 | Z) = logistic(alpha_j + beta_j * Z)

# Set up parameters for logistic regression (one for each treatment)
alpha <- runif(p, min = -1, max = 1)  # intercepts
beta <- runif(p, min = 0.5, max = 2)  # slopes (positive to ensure confounding)

# Generate binary treatments
X <- matrix(0, n, p)
colnames(X) <- paste0("X", 1:p)

for (j in 1:p) {
  # Linear predictor
  eta <- alpha[j] + beta[j] * Z + rnorm(n, mean = 0, sd = 0.1)  # small noise

  # Probability of treatment
  prob <- 1 / (1 + exp(-eta))

  # Generate binary treatment
  X[, j] <- rbinom(n, size = 1, prob = prob)
}

cat("Treatment generation summary:\n")
cat("Treatment prevalence (proportion of 1s):\n")
print(round(colMeans(X), 3))
cat("\n")

# ============================================================================
# 3. Generate Continuous Response (with confounding)
# ============================================================================

# Set causal effects
# True causal treatments: X5, X10, X15 with different effect sizes
delta <- rep(0, p)
delta[5] <- 2.5    # X5 has positive effect
delta[10] <- -2.0  # X10 has negative effect
delta[15] <- 1.5   # X15 has positive effect

# Intercept and confounder effect
intercept <- 2.0
gamma_Z <- 1.5  # Confounding effect

# Generate response with Gaussian noise
epsilon <- rnorm(n, mean = 0, sd = 1)
y <- intercept + X %*% delta + gamma_Z * Z + epsilon

cat("Response generation summary:\n")
cat(sprintf("Mean(y) = %.3f, SD(y) = %.3f\n", mean(y), sd(y)))
cat(sprintf("True causal treatments: X5 (delta=%.2f), X10 (delta=%.2f), X15 (delta=%.2f)\n",
            delta[5], delta[10], delta[15]))
cat(sprintf("Confounder effect: gamma_Z = %.2f\n\n", gamma_Z))

# ============================================================================
# 4. Compute Inverse Propensity Weights (IPW)
# ============================================================================

# Estimate propensity scores for each treatment using logistic regression
W <- matrix(1, n, p)

for (j in 1:p) {
  # Fit logistic regression: X_j ~ Z
  glm_fit <- glm(X[, j] ~ Z, family = binomial(link = "logit"))

  # Predicted propensity scores
  ps <- predict(glm_fit, type = "response")

  # Compute IPW: w_ij = 1/ps_i when X_ij = 1, and 1/(1-ps_i) when X_ij = 0
  W[, j] <- ifelse(X[, j] == 1, 1 / ps, 1 / (1 - ps))

  # Truncate extreme weights to avoid numerical issues
  W[, j] <- pmin(W[, j], quantile(W[, j], 0.99))
  W[, j] <- pmax(W[, j], quantile(W[, j], 0.01))
}

cat("IPW summary:\n")
cat("Weight ranges by treatment (showing first 10):\n")
for (j in 1:min(10, p)) {
  cat(sprintf("  X%d: [%.3f, %.3f], mean = %.3f\n",
              j, min(W[, j]), max(W[, j]), mean(W[, j])))
}
cat("\n")

# ============================================================================
# 5. Source Required Functions
# ============================================================================

# Assuming we're in the examples directory
source("../R/clamp_binary.R")
source("../R/clamp_update_each_effect_binary.R")
source("../R/ipw_single_effect_regression_binary.R")
source("../R/bootstrap_ipw_variance_binary.R")
source("../R/estimate_average_treatment_effect_binary.R")
source("../R/compute_colstats.R")
source("../R/initialize.R")
source("../R/elbo.R")
source("../R/robust_importance_weights.R")
source("../R/remove_abnormal_subjects.R")
source("../R/check_abnormal_subjects.R")
source("../R/sparse_multiplication.R")
source("../R/estimate_residual_variance.R")
source("../R/optimize_prior_variance.R")
source("../R/model_weighted_linear.R")
source("../R/clamp_utils.R")
source("../R/clamp_summarize_coefficients.R")
source("../R/summary.clamp.R")

# Fit clamp_binary model
fit <- clamp_binary(
  X = X,
  y = y,
  W = W,
  maxL = 10,  # Increased number of effects for more variables
  intercept = TRUE,
  mle_estimator = "mHT",
  verbose = TRUE
)

# ============================================================================
# 6. Display Results
# ============================================================================

cat("\n")
cat("RESULTS:\n")
cat("======================================================================\n")
cat("\n")

# Extract PIPs
pip_values <- fit$pip
names(pip_values) <- colnames(X)

# Sort by PIP
pip_sorted <- sort(pip_values, decreasing = TRUE)

# Mark true causal variables
true_causal <- c(5, 10, 15)

cat("Posterior Inclusion Probabilities (PIP):\n")
for (i in 1:length(pip_sorted)) {
  var_idx <- as.numeric(sub("X", "", names(pip_sorted)[i]))
  marker <- if (var_idx %in% true_causal) " ***" else ""
  cat(sprintf("  %s: %.4f%s\n", names(pip_sorted)[i], pip_sorted[i], marker))
}
cat("  (*** indicates true causal variable)\n\n")

# Posterior mean effects
posterior_mean <- colSums(fit$alpha * fit$mu)
names(posterior_mean) <- colnames(X)
pm_sorted <- posterior_mean[names(pip_sorted)]

cat("Posterior Mean Causal Effects:\n")
for (i in 1:length(pm_sorted)) {
  var_idx <- as.numeric(sub("X", "", names(pm_sorted)[i]))
  marker <- if (var_idx %in% true_causal) " ***" else ""
  true_val <- delta[var_idx]
  cat(sprintf("  %s: %.4f (true: %.2f)%s\n",
              names(pm_sorted)[i], pm_sorted[i], true_val, marker))
}
cat("\n")

# Credible sets
if (!is.null(fit$sets$cs)) {
  cat("Credible Sets:\n")
  for (i in 1:length(fit$sets$cs)) {
    cs_vars <- paste0("X", fit$sets$cs[[i]])
    cat(sprintf("  CS%d: {%s}\n", i, paste(cs_vars, collapse = ", ")))
  }
  cat("\n")
}

# Model statistics
cat("Model Fit Statistics:\n")
cat(sprintf("  Number of iterations: %d\n", fit$niter))
cat(sprintf("  Converged: %s\n", fit$converged))
cat(sprintf("  Final ELBO: %.4f\n", fit$elbo[length(fit$elbo)]))
cat(sprintf("  Estimated residual variance: %.4f (true: %.2f)\n",
            fit$sigma2, 1.0))
cat(sprintf("  Intercept: %.4f (true: %.2f)\n", fit$intercept, intercept))
cat("\n")

# Detection summary
detected <- sum(pip_values > 0.5)
true_detected <- sum(pip_values[true_causal] > 0.5)
false_positives <- sum(pip_values[-true_causal] > 0.5)

cat("Detection Summary:\n")
cat(sprintf("  True causal variables: %d\n", length(true_causal)))
cat(sprintf("  Correctly detected (PIP > 0.5): %d\n", true_detected))
cat(sprintf("  False positives (PIP > 0.5): %d\n", false_positives))
cat("\n")

cat("======================================================================\n")
cat("Test completed successfully!\n")

# Display summary table
clamp_summarize_coefficients(fit, top_n = 10)
