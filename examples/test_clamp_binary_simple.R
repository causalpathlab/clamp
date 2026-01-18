#!/usr/bin/env Rscript
# Simple test example for clamp_binary()

# Clean workspace
rm(list = ls())

# Load required packages
library(Matrix)
library(matrixStats)
library(sandwich)

# Source all required functions
cat("Loading functions...\n")
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
source("../R/summary.clamp.R")

set.seed(123)

# ============================================================================
# Simulate Data
# ============================================================================
cat("\n=== Simulating Data ===\n")

n <- 1000  # samples
p <- 10    # treatments

# Confounder
Z <- rnorm(n)
cat(sprintf("Generated confounder Z (n=%d)\n", n))

# Binary treatments confounded by Z
# P(X[,j] = 1 | Z) = logistic(alpha_j + beta_j * Z)
X <- matrix(0, n, p)
colnames(X) <- paste0("X", 1:p)

for (j in 1:p) {
  alpha_j <- runif(1, -1, 1)
  beta_j <- runif(1, 0.5, 2)
  prob <- plogis(alpha_j + beta_j * Z)
  X[, j] <- rbinom(n, 1, prob)
}
cat(sprintf("Generated %d binary treatments\n", p))
cat("Treatment prevalence:", paste(sprintf("%.2f", colMeans(X)), collapse=", "), "\n")

# Response: y = 2 + 1.5*Z + 3*X3 - 2.5*X7 + epsilon
# Only X3 and X7 are truly causal
y <- 2 + 1.5 * Z + 3 * X[, 3] - 2.5 * X[, 7] + rnorm(n, sd = 1)
cat(sprintf("Generated response y (mean=%.2f, sd=%.2f)\n", mean(y), sd(y)))
cat("TRUE CAUSAL VARIABLES: X3 (effect=+3.0), X7 (effect=-2.5)\n")

# Inverse Probability Weights
PS_raw <- matrix(0, n, p)  # Store raw propensity scores
colnames(PS_raw) <- paste0("X", 1:p)

for (j in 1:p) {
  glm_fit <- glm(X[, j] ~ Z, family = binomial())
  ps <- fitted(glm_fit)
  PS_raw[, j] <- ps
}

# Truncate propensity scores to [0.1, 0.9]
PS <- pmax(0.1, pmin(0.9, PS_raw))
colnames(PS) <- paste0("X", 1:p)

# Construct weight matrix from truncated propensity scores
W <- ifelse(X == 1, 1 / PS, 1 / (1 - PS))
colnames(W) <- paste0("X", 1:p)

cat("Computed propensity scores and IPW matrix\n")
cat(sprintf("Truncated propensity scores to [0.1, 0.9]\n"))

# ============================================================================
# Check Propensity Score Distributions
# ============================================================================
cat("\n=== Propensity Score Analysis ===\n")

# Create histograms for raw propensity scores
pdf("propensity_score_raw.pdf", width = 12, height = 8)
par(mfrow = c(2, 5))
for (j in 1:p) {
  hist(PS_raw[, j],
       main = paste0("Raw PS for ", colnames(PS_raw)[j]),
       xlab = "Propensity Score",
       col = "#71bdf4",
       breaks = 30,
       xlim = c(0, 1))
  abline(v = c(0.1, 0.9), col = "red", lwd = 2, lty = 2)
}
dev.off()
cat("Saved raw propensity score histograms to 'propensity_score_raw.pdf'\n")

# Create histograms for truncated propensity scores
pdf("propensity_score_truncated.pdf", width = 12, height = 8)
par(mfrow = c(2, 5))
for (j in 1:p) {
  hist(PS[, j],
       main = paste0("Truncated PS for ", colnames(PS)[j]),
       xlab = "Propensity Score",
       col = "#90EE90",
       breaks = 30,
       xlim = c(0, 1))
  abline(v = c(0.1, 0.9), col = "red", lwd = 2, lty = 2)
}
dev.off()
cat("Saved truncated propensity score histograms to",
    "'propensity_score_truncated.pdf'\n")

# Create histograms for weights
pdf("weight_distributions.pdf", width = 12, height = 8)
par(mfrow = c(2, 5))
for (j in 1:p) {
  hist(W[, j],
       main = paste0("Weights for ", colnames(W)[j]),
       xlab = "Weight",
       col = "#FFB6C1",
       breaks = 30)
}
dev.off()
cat("Saved weight histograms to 'weight_distributions.pdf'\n")

# Report extreme propensity scores in raw data
cat("\n=== Extreme Propensity Score Summary (Raw) ===\n")
cat("(Extreme: PS < 0.1 or PS > 0.9 before truncation)\n\n")

total_extreme <- 0
for (j in 1:p) {
  n_low <- sum(PS_raw[, j] < 0.1)
  n_high <- sum(PS_raw[, j] > 0.9)
  n_extreme <- n_low + n_high
  total_extreme <- total_extreme + n_extreme

  pct_extreme <- 100 * n_extreme / n

  cat(sprintf("%s: %d extreme (%.1f%%) [%d < 0.1, %d > 0.9]\n",
              colnames(PS_raw)[j], n_extreme, pct_extreme, n_low, n_high))
}

cat(sprintf("\nTotal extreme values: %d out of %d (%.1f%%)\n",
            total_extreme, n * p, 100 * total_extreme / (n * p)))

# Summary statistics for raw propensity scores
cat("\n=== Raw Propensity Score Summary Statistics ===\n")
ps_raw_summary <- data.frame(
  Variable = colnames(PS_raw),
  Min = apply(PS_raw, 2, min),
  Median = apply(PS_raw, 2, median),
  Mean = apply(PS_raw, 2, mean),
  Max = apply(PS_raw, 2, max),
  SD = apply(PS_raw, 2, sd)
)
print(ps_raw_summary, digits = 3, row.names = FALSE)

# Summary statistics for truncated propensity scores
cat("\n=== Truncated Propensity Score Summary Statistics ===\n")
ps_summary <- data.frame(
  Variable = colnames(PS),
  Min = apply(PS, 2, min),
  Median = apply(PS, 2, median),
  Mean = apply(PS, 2, mean),
  Max = apply(PS, 2, max),
  SD = apply(PS, 2, sd)
)
print(ps_summary, digits = 3, row.names = FALSE)

# Summary statistics for weights
cat("\n=== Weight Summary Statistics ===\n")
weight_summary <- data.frame(
  Variable = colnames(W),
  Min = apply(W, 2, min),
  Median = apply(W, 2, median),
  Mean = apply(W, 2, mean),
  Max = apply(W, 2, max),
  SD = apply(W, 2, sd)
)
print(weight_summary, digits = 3, row.names = FALSE)



# ============================================================================
# Run clamp_binary
# ============================================================================
cat("\n=== Running clamp_binary ===\n")

fit <- clamp_binary(
  X = X,
  y = y,
  W = W,
  maxL = 5,
  mle_estimator = "mHT",
  mle_variance_estimator = "bootstrap",
  nboots = 50,
  standardize = FALSE,
  intercept = TRUE,
  estimate_residual_variance = TRUE,
  estimate_prior_variance = TRUE,
  max_iter = 100,
  tol = 1e-2,
  verbose = FALSE
)

# ============================================================================
# Results
# ============================================================================
cat("\n=== RESULTS ===\n")

# PIPs
cat("\nPosterior Inclusion Probabilities:\n")
pip_df <- data.frame(
  Variable = names(fit$pip),
  PIP = fit$pip,
  True = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE)
)
pip_df <- pip_df[order(-pip_df$PIP), ]
rownames(pip_df) <- NULL
print(pip_df, digits = 4)

# Effects
cat("\nPosterior Mean Effects:\n")
effects <- colSums(fit$alpha * fit$mu)
true_effects <- c(0, 0, 3, 0, 0, 0, -2.5, 0, 0, 0)
names(true_effects) <- paste0("X", 1:10)
effect_df <- data.frame(
  Variable = names(effects),
  Estimated = effects,
  True = true_effects
)
effect_df <- effect_df[order(-abs(effect_df$Estimated)), ]
rownames(effect_df) <- NULL
print(effect_df, digits = 3)

# Summary
cat("\nModel Summary:\n")
cat(sprintf("  Iterations: %d (converged: %s)\n", fit$niter, fit$converged))
cat(sprintf("  Residual variance: %.3f (true: 1.0)\n", fit$sigma2))
cat(sprintf("  Intercept: %.3f (true: 2.0)\n", fit$intercept))

detected <- names(fit$pip)[fit$pip > 0.5]
cat(sprintf("\nDetected variables (PIP > 0.5): %s\n",
            ifelse(length(detected) > 0, paste(detected, collapse = ", "), "None")))

# Check success
if ("X3" %in% detected && "X7" %in% detected) {
  cat("\n✓ SUCCESS: Both true causal variables detected!\n")
} else {
  cat("\n✗ WARNING: Not all true variables detected.\n")
}
