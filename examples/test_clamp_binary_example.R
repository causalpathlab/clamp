#!/usr/bin/env Rscript
# Test example for clamp_binary() function
# Simulates confounded binary treatments and tests causal variable selection

# Load required packages
library(Matrix)
library(matrixStats)
library(sandwich)

set.seed(123)

# ============================================================================
# 1. Simulate Data
# ============================================================================

n <- 1000  # number of samples
p <- 10    # number of treatments

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
X <- matrix(0, nrow = n, ncol = p)
colnames(X) <- paste0("X", 1:p)

for (j in 1:p) {
  # Probability of treatment
  prob_treat <- plogis(alpha[j] + beta[j] * Z)
  # Generate binary treatment
  X[, j] <- rbinom(n, size = 1, prob = prob_treat)
}

cat("Treatment generation summary:\n")
cat("Treatment prevalence (proportion of 1s):\n")
print(round(colMeans(X), 3))
cat("\n")

# ============================================================================
# 3. Generate Response (only 2 true causal treatments + confounder)
# ============================================================================
# True causal treatments: X3 and X7
# y = gamma_0 + gamma_Z * Z + delta_3 * X3 + delta_7 * X7 + epsilon

true_causal_vars <- c(3, 7)  # X3 and X7 are truly causal
gamma_0 <- 2.0               # intercept
gamma_Z <- 1.5               # confounder effect on outcome
delta_3 <- 3.0               # true causal effect of X3
delta_7 <- -2.5              # true causal effect of X7
sigma <- 1.0                 # residual standard deviation

# Generate response
y <- gamma_0 +
     gamma_Z * Z +
     delta_3 * X[, 3] +
     delta_7 * X[, 7] +
     rnorm(n, mean = 0, sd = sigma)

cat("Response generation summary:\n")
cat(sprintf("Mean(y) = %.3f, SD(y) = %.3f\n", mean(y), sd(y)))
cat(sprintf("True causal treatments: X%d (delta=%.2f), X%d (delta=%.2f)\n",
            true_causal_vars[1], delta_3, true_causal_vars[2], delta_7))
cat(sprintf("Confounder effect: gamma_Z = %.2f\n", gamma_Z))
cat("\n")

# ============================================================================
# 4. Compute Inverse Probability Weights (IPW)
# ============================================================================
# For each treatment, estimate propensity score using logistic regression
# Then compute inverse probability weights

W <- matrix(0, nrow = n, ncol = p)
colnames(W) <- paste0("X", 1:p)

for (j in 1:p) {
  # Fit logistic regression: P(X[,j] = 1 | Z)
  glm_fit <- glm(X[, j] ~ Z, family = binomial(link = "logit"))

  # Predicted propensity scores
  ps <- predict(glm_fit, type = "response")

  # Compute IPW: w_ij = X_ij / ps_i + (1 - X_ij) / (1 - ps_i)
  # But in clamp, we use: w_ij = 1/ps_i when X_ij = 1, and 1/(1-ps_i) when X_ij = 0
  W[, j] <- ifelse(X[, j] == 1, 1 / ps, 1 / (1 - ps))

  # Truncate extreme weights to avoid numerical issues
  W[, j] <- pmin(W[, j], quantile(W[, j], 0.99))
  W[, j] <- pmax(W[, j], quantile(W[, j], 0.01))
}

cat("IPW summary:\n")
cat("Weight ranges by treatment:\n")
for (j in 1:p) {
  cat(sprintf("  X%d: [%.3f, %.3f], mean = %.3f\n",
              j, min(W[, j]), max(W[, j]), mean(W[, j])))
}
cat("\n")

# ============================================================================
# 5. Run clamp_binary()
# ============================================================================
cat("Running clamp_binary()...\n")
cat(rep("=", 70), "\n", sep = "")

# Source the clamp_binary function and dependencies
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

# Fit clamp_binary model
fit <- clamp_binary(
  X = X,
  y = y,
  W = W,
  maxL = 5,                           # Allow up to 5 effects
  mle_estimator = "mHT",              # modified Horvitz-Thompson estimator
  mle_variance_estimator = "bootstrap", # bootstrap variance estimation
  nboots = 50,                        # reduce for faster computation
  standardize = FALSE,                # Don't standardize for binary data
  intercept = TRUE,
  estimate_residual_variance = TRUE,
  estimate_prior_variance = TRUE,
  estimate_prior_method = "optim",
  max_iter = 100,
  tol = 1e-2,
  verbose = TRUE
)

cat(rep("=", 70), "\n", sep = "")
cat("\n")

# ============================================================================
# 6. Display Results
# ============================================================================
cat("RESULTS:\n")
cat(rep("=", 70), "\n", sep = "")

# Posterior Inclusion Probabilities (PIP)
cat("\nPosterior Inclusion Probabilities (PIP):\n")
pip_sorted <- sort(fit$pip, decreasing = TRUE)
for (i in 1:p) {
  var_name <- names(pip_sorted)[i]
  pip_val <- pip_sorted[i]
  is_true <- var_name %in% paste0("X", true_causal_vars)
  marker <- ifelse(is_true, " ***", "")
  cat(sprintf("  %s: %.4f%s\n", var_name, pip_val, marker))
}
cat("  (*** indicates true causal variable)\n\n")

# Posterior mean effects
cat("Posterior Mean Causal Effects:\n")
posterior_effects <- colSums(fit$alpha * fit$mu)
effects_sorted <- sort(abs(posterior_effects), decreasing = TRUE)
for (i in 1:p) {
  var_name <- names(effects_sorted)[i]
  effect_val <- posterior_effects[var_name]
  is_true <- var_name %in% paste0("X", true_causal_vars)

  if (var_name == "X3") {
    true_val <- delta_3
  } else if (var_name == "X7") {
    true_val <- delta_7
  } else {
    true_val <- 0.0
  }

  marker <- ifelse(is_true, " ***", "")
  cat(sprintf("  %s: %.4f (true: %.2f)%s\n", var_name, effect_val, true_val, marker))
}
cat("\n")

# Credible Sets
if (!is.null(fit$sets) && length(fit$sets$cs) > 0) {
  cat("Credible Sets:\n")
  for (i in 1:length(fit$sets$cs)) {
    cs_vars <- colnames(X)[fit$sets$cs[[i]]]
    cat(sprintf("  CS%d: {%s}\n", i, paste(cs_vars, collapse = ", ")))
  }
  cat("\n")
} else {
  cat("No credible sets found.\n\n")
}

# Model fit statistics
cat("Model Fit Statistics:\n")
cat(sprintf("  Number of iterations: %d\n", fit$niter))
cat(sprintf("  Converged: %s\n", fit$converged))
cat(sprintf("  Final ELBO: %.4f\n", tail(fit$elbo, 1)))
cat(sprintf("  Estimated residual variance: %.4f (true: %.2f)\n",
            fit$sigma2, sigma^2))
cat(sprintf("  Intercept: %.4f (true: %.2f)\n", fit$intercept, gamma_0))
cat("\n")

# Check if true variables were detected
detected_vars <- names(fit$pip)[fit$pip > 0.5]
true_var_names <- paste0("X", true_causal_vars)
correctly_detected <- sum(true_var_names %in% detected_vars)
false_positives <- sum(!(detected_vars %in% true_var_names))

cat("Detection Summary:\n")
cat(sprintf("  True causal variables: %d\n", length(true_causal_vars)))
cat(sprintf("  Correctly detected (PIP > 0.5): %d\n", correctly_detected))
cat(sprintf("  False positives (PIP > 0.5): %d\n", false_positives))
cat("\n")

cat(rep("=", 70), "\n", sep = "")
cat("Test completed successfully!\n")
