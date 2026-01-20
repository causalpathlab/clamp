#!/usr/bin/env Rscript
# Test example for clamp_binary() with IMBALANCED binary treatments
# Using OVERLAP WEIGHTS instead of IPW
# Some treatments have extreme imbalance (80:20 or 90:10)

# Load required packages
library(Matrix)
library(matrixStats)
library(sandwich)

set.seed(456)

# ============================================================================
# 1. Simulate Data
# ============================================================================

n <- 1000  # number of samples
p <- 20    # number of treatments

# Generate 1-dimensional confounder
Z <- rnorm(n, mean = 0, sd = 1)

# ============================================================================
# 2. Generate Binary Treatments with IMBALANCED distributions
# ============================================================================
# Create different levels of imbalance:
# - Variables 1-5: Moderately imbalanced (~30% treatment)
# - Variables 6-10: Imbalanced (~20% treatment)
# - Variables 11-15: Highly imbalanced (~10% treatment)
# - Variables 16-20: Extremely imbalanced (~5% treatment)

X <- matrix(0, n, p)
colnames(X) <- paste0("X", 1:p)

# Define intercepts to control treatment prevalence
# Lower intercept = lower probability of treatment
alpha_base <- c(
  rep(-0.5, 5),   # ~30% treatment (X1-X5)
  rep(-1.2, 5),   # ~20% treatment (X6-X10)
  rep(-2.0, 5),   # ~10% treatment (X11-X15)
  rep(-2.8, 5)    # ~5% treatment (X16-X20)
)

# Confounding effects (vary across treatments)
beta <- runif(p, min = 0.3, max = 1.5)

for (j in 1:p) {
  # Linear predictor with confounding
  eta <- alpha_base[j] + beta[j] * Z

  # Probability of treatment
  prob <- 1 / (1 + exp(-eta))

  # Generate binary treatment
  X[, j] <- rbinom(n, size = 1, prob = prob)
}

cat("Treatment generation summary:\n")
cat("Treatment prevalence (proportion of 1s):\n")
prevalence <- colMeans(X)
print(round(prevalence, 3))
cat("\n")

# Display prevalence by group
cat("Prevalence by design group:\n")
cat(sprintf("  X1-X5 (target ~30%%): mean = %.1f%%\n",
            mean(prevalence[1:5]) * 100))
cat(sprintf("  X6-X10 (target ~20%%): mean = %.1f%%\n",
            mean(prevalence[6:10]) * 100))
cat(sprintf("  X11-X15 (target ~10%%): mean = %.1f%%\n",
            mean(prevalence[11:15]) * 100))
cat(sprintf("  X16-X20 (target ~5%%): mean = %.1f%%\n",
            mean(prevalence[16:20]) * 100))
cat("\n")

# ============================================================================
# 3. Generate Response with TRUE causal effects
# ============================================================================
# True causal treatments: X3 (moderately imbalanced) and X12 (highly imbalanced)

true_causal_vars <- c(3, 12)
gamma_0 <- 2.0
gamma_Z <- 1.5
delta_3 <- 3.0    # Effect from moderately imbalanced treatment
delta_12 <- -2.5  # Effect from highly imbalanced treatment
sigma <- 1.0
eps <- rnorm(n, mean = 0, sd = sigma)

# Create delta vector for all variables
delta <- rep(0, p)
delta[3] <- delta_3
delta[12] <- delta_12

# Generate response
y <- gamma_0 + gamma_Z * Z + delta_3 * X[, 3] + delta_12 * X[, 12] + eps

cat("Response generation summary:\n")
cat(sprintf("Mean(y) = %.3f, SD(y) = %.3f\n", mean(y), sd(y)))
cat(sprintf("True causal treatments: X3 (delta=%.2f, prevalence=%.1f%%), X12 (delta=%.2f, prevalence=%.1f%%)\n",
            delta_3, prevalence[3]*100, delta_12, prevalence[12]*100))
cat(sprintf("Confounder effect: gamma_Z = %.2f\n\n", gamma_Z))

# ============================================================================
# 4. Compute OVERLAP WEIGHTS
# ============================================================================

W <- matrix(0, nrow = n, ncol = p)
colnames(W) <- paste0("X", 1:p)

# Matrix to store propensity scores
PS <- matrix(0, nrow = n, ncol = p)
colnames(PS) <- paste0("X", 1:p)

cat("Computing OVERLAP WEIGHTS...\n")
cat("Overlap weight formula: w = ps * (1 - ps)\n\n")

for (j in 1:p) {
  # Fit logistic regression: P(X[,j] = 1 | Z)
  glm_fit <- glm(X[, j] ~ Z, family = binomial(link = "logit"))

  # Predicted propensity scores
  ps <- predict(glm_fit, type = "response")
  PS[, j] <- ps

  # Compute OVERLAP WEIGHTS: w = ps * (1 - ps)
  W[, j] <- ps * (1 - ps)
}

# ============================================================================
# 4.5 Display Propensity Score Distributions
# ============================================================================

cat("PROPENSITY SCORE DISTRIBUTIONS\n")
cat(rep("=", 70), "\n", sep = "")
cat("\n")

for (j in 1:p) {
  ps <- PS[, j]
  cat(sprintf("X%d (prevalence = %.1f%%):\n", j, prevalence[j]*100))
  cat(sprintf("  Propensity scores: min=%.4f, Q1=%.4f, median=%.4f, Q3=%.4f, max=%.4f\n",
              min(ps), quantile(ps, 0.25), median(ps),
              quantile(ps, 0.75), max(ps)))
  cat(sprintf("  Overlap weights: min=%.4f, mean=%.4f, max=%.4f\n",
              min(W[, j]), mean(W[, j]), max(W[, j])))
  cat("\n")
}

cat("\nOverlap Weight Summary (by imbalance group):\n")
for (group in list(1:5, 6:10, 11:15, 16:20)) {
  group_name <- sprintf("X%d-X%d", min(group), max(group))
  group_weights <- W[, group]
  cat(sprintf("  %s: mean weight = %.4f (SD = %.4f)\n",
              group_name, mean(group_weights), sd(as.vector(group_weights))))
}
cat("\n")

# ============================================================================
# 5. Run clamp_binary()
# ============================================================================

cat("Running clamp_binary() with OVERLAP WEIGHTS...\n")
cat(rep("=", 70), "\n", sep = "")

# Source required functions
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
  maxL = 5,
  mle_estimator = "mHT",
  mle_variance_estimator = "bootstrap",
  nboots = 50,
  intercept = TRUE,
  estimate_residual_variance = TRUE,
  estimate_prior_variance = TRUE,
  max_iter = 100,
  tol = 1e-2,
  verbose = TRUE
)

cat(rep("=", 70), "\n", sep = "")

# ============================================================================
# 6. Display Results
# ============================================================================

cat("\n")
cat("RESULTS:\n")
cat(rep("=", 70), "\n", sep = "")
cat("\n")

# Extract PIPs
pip_values <- fit$pip
names(pip_values) <- colnames(X)
pip_sorted <- sort(pip_values, decreasing = TRUE)

cat("Posterior Inclusion Probabilities (PIP):\n")
for (i in 1:length(pip_sorted)) {
  var_idx <- as.numeric(sub("X", "", names(pip_sorted)[i]))
  marker <- if (var_idx %in% true_causal_vars) " ***" else ""
  prev_marker <- sprintf("(prev=%.1f%%)", prevalence[var_idx]*100)
  cat(sprintf("  %s: %.4f %s %s\n",
              names(pip_sorted)[i], pip_sorted[i], prev_marker, marker))
}
cat("  (*** indicates true causal variable)\n\n")

# Posterior mean effects
posterior_mean <- colSums(fit$alpha * fit$mu)
names(posterior_mean) <- colnames(X)
pm_sorted <- posterior_mean[names(pip_sorted)]

cat("Posterior Mean Causal Effects:\n")
for (i in 1:length(pm_sorted)) {
  var_idx <- as.numeric(sub("X", "", names(pm_sorted)[i]))
  marker <- if (var_idx %in% true_causal_vars) " ***" else ""
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
            fit$sigma2, sigma^2))
cat(sprintf("  Intercept: %.4f (true: %.2f)\n", fit$intercept, gamma_0))
cat("\n")

# Detection summary
detected <- sum(pip_values > 0.5)
true_detected <- sum(pip_values[true_causal_vars] > 0.5)
false_positives <- sum(pip_values[-true_causal_vars] > 0.5)

cat("Detection Summary:\n")
cat(sprintf("  True causal variables: %d\n", length(true_causal_vars)))
cat(sprintf("  Correctly detected (PIP > 0.5): %d\n", true_detected))
cat(sprintf("  False positives (PIP > 0.5): %d\n", false_positives))
cat("\n")

cat(rep("=", 70), "\n", sep = "")
cat("Test completed successfully!\n\n")

# Display summary table
cat("clamp_binary Results:\n")
clamp_summarize_coefficients(fit, top_n = 15)

# ============================================================================
# 7. Analysis of Imbalanced Variables
# ============================================================================

cat("\n\n")
cat(rep("=", 70), "\n", sep = "")
cat("ANALYSIS: Effect of Treatment Imbalance on Detection\n")
cat(rep("=", 70), "\n", sep = "")
cat("\n")

# Create analysis table
imbalance_analysis <- data.frame(
  Variable = colnames(X),
  Prevalence = prevalence * 100,
  Mean_PS = colMeans(PS),
  Mean_Weight = colMeans(W),
  PIP = fit$pip,
  Effect = colSums(fit$alpha * fit$mu),
  True_Effect = delta,
  Is_True = ifelse(1:p %in% true_causal_vars, "YES", "NO")
)

# Sort by prevalence
imbalance_analysis <- imbalance_analysis[order(imbalance_analysis$Prevalence,
                                                decreasing = TRUE), ]

cat("Effect of treatment prevalence on detection:\n\n")
cat("Variables sorted by prevalence (highest to lowest):\n\n")

# Create a copy with rounded numeric columns for display
imbalance_display <- imbalance_analysis
numeric_cols <- sapply(imbalance_display, is.numeric)
imbalance_display[numeric_cols] <- lapply(imbalance_display[numeric_cols], round, 4)
print(imbalance_display, row.names = FALSE)

cat("\n\nKey Observations:\n")
cat(sprintf("  - X3 (%.1f%% prevalence): PIP = %.4f, Effect = %.4f (true: %.2f)\n",
            prevalence[3]*100, fit$pip[3],
            colSums(fit$alpha * fit$mu)[3], delta[3]))
cat(sprintf("  - X12 (%.1f%% prevalence, highly imbalanced): PIP = %.4f, Effect = %.4f (true: %.2f)\n",
            prevalence[12]*100, fit$pip[12],
            colSums(fit$alpha * fit$mu)[12], delta[12]))
cat("\n")

cat(rep("=", 70), "\n", sep = "")
cat("ANALYSIS COMPLETE\n")
cat(rep("=", 70), "\n", sep = "")

# ============================================================================
# 8. Visualize Propensity Score Distributions
# ============================================================================

# Create plots subfolder
plots_dir <- "imbalanced-study-plots"
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
  cat(sprintf("\nCreated directory: %s\n", plots_dir))
} else {
  cat(sprintf("\nUsing existing directory: %s\n", plots_dir))
}

cat("\nCreating histograms for propensity scores...\n")

# Save histograms to PDF
pdf(file.path(plots_dir, "propensity_score_histograms.pdf"),
    width = 15, height = 12)

# Set up plotting area: 4x5 grid for 20 variables
par(mfrow = c(4, 5), mar = c(3, 3, 2, 1), oma = c(0, 0, 2, 0))

for (j in 1:p) {
  hist(PS[, j],
       main = sprintf("X%d (prev=%.1f%%)", j, prevalence[j]*100),
       xlab = "Propensity Score",
       ylab = "Frequency",
       col = "lightblue",
       border = "darkblue",
       breaks = 30,
       xlim = c(0, 1))

  # Add vertical line at mean
  abline(v = mean(PS[, j]), col = "red", lwd = 2, lty = 2)

  # Highlight true causal variables with different color
  if (j %in% true_causal_vars) {
    box(col = "red", lwd = 3)
  }
}

# Add overall title
mtext("Propensity Score Distributions for All Treatment Variables",
      outer = TRUE, cex = 1.2, font = 2)

# Reset plotting parameters
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, oma = c(0, 0, 0, 0))
dev.off()

cat(sprintf("Histograms saved to: %s\n",
            file.path(plots_dir, "propensity_score_histograms.pdf")))
cat("Red boxes indicate true causal variables (X3 and X12)\n")
cat("Red dashed lines show mean propensity scores\n")

# ============================================================================
# 9. Display ALL variables using clamp_summarize_coefficients
# ============================================================================

cat("\n\n")
cat(rep("=", 70), "\n", sep = "")
cat("CLAMP Results for ALL Variables\n")
cat(rep("=", 70), "\n", sep = "")
cat("\n")

clamp_summarize_coefficients(fit, top_n = p)

# ============================================================================
# 10. Comparison with susieR::susie()
# ============================================================================

cat("\n\n")
cat(rep("=", 70), "\n", sep = "")
cat("COMPARISON WITH susieR::susie()\n")
cat(rep("=", 70), "\n", sep = "")
cat("\n")

library(susieR)

# Run standard SuSiE (without IPW weights)
cat("Running susieR::susie() (standard method, no weights)...\n\n")
fit_susie <- susie(X, y, L = 5, intercept = TRUE,
                   estimate_residual_variance = TRUE,
                   estimate_prior_variance = TRUE,
                   max_iter = 100, tol = 1e-2, verbose = FALSE)

cat("SuSiE Results for ALL Variables:\n")
cat(rep("-", 70), "\n", sep = "")

# Extract results from susieR
susie_pip <- fit_susie$pip
susie_effect <- coef(fit_susie)[-1]  # Remove intercept

# Compute posterior standard deviation for each variable
susie_sd <- sqrt(colSums(fit_susie$alpha * fit_susie$mu2) -
                 colSums(fit_susie$alpha * fit_susie$mu)^2)

# Create summary table for SuSiE
susie_summary <- data.frame(
  variable = colnames(X),
  pip = susie_pip,
  posterior_mean = susie_effect,
  posterior_sd = susie_sd
)

# Sort by PIP
susie_summary <- susie_summary[order(susie_summary$pip, decreasing = TRUE), ]
print(susie_summary, row.names = FALSE, digits = 4)

cat("\n")
cat(rep("-", 70), "\n", sep = "")
cat("Detection Performance Comparison (PIP > 0.5):\n")
cat(rep("-", 70), "\n", sep = "")
cat("\n")

clamp_pip <- fit$pip
clamp_detected <- sum(clamp_pip > 0.5)
clamp_true_detected <- sum(clamp_pip[true_causal_vars] > 0.5)
clamp_false_pos <- sum(clamp_pip[-true_causal_vars] > 0.5)

susie_detected <- sum(susie_pip > 0.5)
susie_true_detected <- sum(susie_pip[true_causal_vars] > 0.5)
susie_false_pos <- sum(susie_pip[-true_causal_vars] > 0.5)

cat("CLAMP (with overlap weights):\n")
cat(sprintf("  Total detected: %d\n", clamp_detected))
cat(sprintf("  True positives: %d/%d\n", clamp_true_detected,
            length(true_causal_vars)))
cat(sprintf("  False positives: %d\n", clamp_false_pos))
cat("\n")

cat("SuSiE (standard, no weights):\n")
cat(sprintf("  Total detected: %d\n", susie_detected))
cat(sprintf("  True positives: %d/%d\n", susie_true_detected,
            length(true_causal_vars)))
cat(sprintf("  False positives: %d\n", susie_false_pos))
cat("\n")

# Effect estimation comparison
clamp_effect <- colSums(fit$alpha * fit$mu)
cat(rep("-", 70), "\n", sep = "")
cat("Effect Estimation for True Causal Variables:\n")
cat(rep("-", 70), "\n", sep = "")
cat("\n")

for (idx in true_causal_vars) {
  cat(sprintf("X%d (true = %.2f, prevalence = %.1f%%):\n",
              idx, delta[idx], prevalence[idx] * 100))
  cat(sprintf("  CLAMP: %.4f (error: %+.4f, PIP: %.4f)\n",
              clamp_effect[idx], clamp_effect[idx] - delta[idx],
              clamp_pip[idx]))
  cat(sprintf("  SuSiE: %.4f (error: %+.4f, PIP: %.4f)\n",
              susie_effect[idx], susie_effect[idx] - delta[idx],
              susie_pip[idx]))
  cat("\n")
}

# Credible sets
cat(rep("-", 70), "\n", sep = "")
cat("Credible Sets:\n")
cat(rep("-", 70), "\n", sep = "")
cat("\n")

cat("CLAMP:\n")
if (!is.null(fit$sets$cs)) {
  for (i in 1:length(fit$sets$cs)) {
    cs_vars <- paste0("X", fit$sets$cs[[i]])
    cat(sprintf("  CS%d: {%s}\n", i, paste(cs_vars, collapse = ", ")))
  }
} else {
  cat("  No credible sets\n")
}
cat("\n")

cat("SuSiE:\n")
if (!is.null(fit_susie$sets$cs)) {
  for (i in 1:length(fit_susie$sets$cs)) {
    cs_vars <- paste0("X", fit_susie$sets$cs[[i]])
    cat(sprintf("  CS%d: {%s}\n", i, paste(cs_vars, collapse = ", ")))
  }
} else {
  cat("  No credible sets\n")
}
cat("\n")

cat(rep("=", 70), "\n", sep = "")

# ============================================================================
# 11. Comparison with LASSO (glmnet)
# ============================================================================

cat("\n\n")
cat(rep("=", 70), "\n", sep = "")
cat("COMPARISON WITH LASSO (glmnet)\n")
cat(rep("=", 70), "\n", sep = "")
cat("\n")

library(glmnet)

# Fit LASSO with cross-validation to select lambda
cat("Running glmnet::cv.glmnet() with alpha=1 (LASSO)...\n")
set.seed(456)  # For reproducibility of CV folds
cv_lasso <- cv.glmnet(X, y, alpha = 1, intercept = TRUE,
                      standardize = TRUE, nfolds = 10)

# Extract coefficients at lambda.min and lambda.1se
lasso_coef_min <- as.vector(coef(cv_lasso, s = "lambda.min"))[-1]
lasso_coef_1se <- as.vector(coef(cv_lasso, s = "lambda.1se"))[-1]
lasso_intercept_min <- as.vector(coef(cv_lasso, s = "lambda.min"))[1]
lasso_intercept_1se <- as.vector(coef(cv_lasso, s = "lambda.1se"))[1]

cat(sprintf("Lambda selected by CV: lambda.min = %.6f, lambda.1se = %.6f\n",
            cv_lasso$lambda.min, cv_lasso$lambda.1se))
cat("\n")

# Create summary tables
cat("LASSO Results (lambda.min - less regularization):\n")
cat(rep("-", 70), "\n", sep = "")

lasso_summary_min <- data.frame(
  variable = colnames(X),
  coefficient = lasso_coef_min,
  selected = ifelse(abs(lasso_coef_min) > 1e-10, "YES", "NO")
)
lasso_summary_min <- lasso_summary_min[order(abs(lasso_summary_min$coefficient),
                                              decreasing = TRUE), ]
print(lasso_summary_min, row.names = FALSE, digits = 4)

cat("\n")
cat("LASSO Results (lambda.1se - more regularization):\n")
cat(rep("-", 70), "\n", sep = "")

lasso_summary_1se <- data.frame(
  variable = colnames(X),
  coefficient = lasso_coef_1se,
  selected = ifelse(abs(lasso_coef_1se) > 1e-10, "YES", "NO")
)
lasso_summary_1se <- lasso_summary_1se[order(abs(lasso_summary_1se$coefficient),
                                              decreasing = TRUE), ]
print(lasso_summary_1se, row.names = FALSE, digits = 4)

cat("\n")
cat(rep("-", 70), "\n", sep = "")
cat("Detection Performance Comparison:\n")
cat(rep("-", 70), "\n", sep = "")
cat("\n")

# Detection statistics
lasso_detected_min <- sum(abs(lasso_coef_min) > 1e-10)
lasso_true_detected_min <- sum(abs(lasso_coef_min[true_causal_vars]) > 1e-10)
lasso_false_pos_min <- lasso_detected_min - lasso_true_detected_min

lasso_detected_1se <- sum(abs(lasso_coef_1se) > 1e-10)
lasso_true_detected_1se <- sum(abs(lasso_coef_1se[true_causal_vars]) > 1e-10)
lasso_false_pos_1se <- lasso_detected_1se - lasso_true_detected_1se

cat("CLAMP (with overlap weights):\n")
cat(sprintf("  Total detected (PIP > 0.5): %d\n", sum(fit$pip > 0.5)))
cat(sprintf("  True positives: %d/%d\n", sum(fit$pip[true_causal_vars] > 0.5),
            length(true_causal_vars)))
cat(sprintf("  False positives: %d\n",
            sum(fit$pip > 0.5) - sum(fit$pip[true_causal_vars] > 0.5)))
cat("\n")

cat("SuSiE (standard, no weights):\n")
cat(sprintf("  Total detected (PIP > 0.5): %d\n", sum(susie_pip > 0.5)))
cat(sprintf("  True positives: %d/%d\n", sum(susie_pip[true_causal_vars] > 0.5),
            length(true_causal_vars)))
cat(sprintf("  False positives: %d\n",
            sum(susie_pip > 0.5) - sum(susie_pip[true_causal_vars] > 0.5)))
cat("\n")

cat("LASSO (lambda.min):\n")
cat(sprintf("  Total detected (non-zero): %d\n", lasso_detected_min))
cat(sprintf("  True positives: %d/%d\n", lasso_true_detected_min,
            length(true_causal_vars)))
cat(sprintf("  False positives: %d\n", lasso_false_pos_min))
cat("\n")

cat("LASSO (lambda.1se):\n")
cat(sprintf("  Total detected (non-zero): %d\n", lasso_detected_1se))
cat(sprintf("  True positives: %d/%d\n", lasso_true_detected_1se,
            length(true_causal_vars)))
cat(sprintf("  False positives: %d\n", lasso_false_pos_1se))
cat("\n")

# Effect estimation comparison
cat(rep("-", 70), "\n", sep = "")
cat("Effect Estimation for True Causal Variables:\n")
cat(rep("-", 70), "\n", sep = "")
cat("\n")

clamp_effect <- colSums(fit$alpha * fit$mu)
for (idx in true_causal_vars) {
  cat(sprintf("X%d (true = %.2f, prevalence = %.1f%%):\n",
              idx, delta[idx], prevalence[idx] * 100))
  cat(sprintf("  CLAMP:           %.4f (error: %+.4f)\n",
              clamp_effect[idx], clamp_effect[idx] - delta[idx]))
  cat(sprintf("  SuSiE:           %.4f (error: %+.4f)\n",
              susie_effect[idx], susie_effect[idx] - delta[idx]))
  cat(sprintf("  LASSO (min):     %.4f (error: %+.4f)\n",
              lasso_coef_min[idx], lasso_coef_min[idx] - delta[idx]))
  cat(sprintf("  LASSO (1se):     %.4f (error: %+.4f)\n",
              lasso_coef_1se[idx], lasso_coef_1se[idx] - delta[idx]))
  cat("\n")
}

# Summary statistics
cat(rep("-", 70), "\n", sep = "")
cat("Model Fit Summary:\n")
cat(rep("-", 70), "\n", sep = "")
cat("\n")

# Compute R-squared for LASSO models
pred_min <- predict(cv_lasso, newx = X, s = "lambda.min")
pred_1se <- predict(cv_lasso, newx = X, s = "lambda.1se")
ss_tot <- sum((y - mean(y))^2)
ss_res_min <- sum((y - pred_min)^2)
ss_res_1se <- sum((y - pred_1se)^2)
r2_min <- 1 - ss_res_min / ss_tot
r2_1se <- 1 - ss_res_1se / ss_tot

cat("CLAMP:\n")
cat(sprintf("  Residual variance: %.4f (true: %.2f)\n",
            fit$sigma2, sigma^2))
cat(sprintf("  Intercept: %.4f (true: %.2f)\n", fit$intercept, gamma_0))
cat("\n")

cat("SuSiE:\n")
cat(sprintf("  Residual variance: %.4f (true: %.2f)\n",
            fit_susie$sigma2, sigma^2))
cat(sprintf("  Intercept: %.4f (true: %.2f)\n", fit_susie$intercept, gamma_0))
cat("\n")

cat("LASSO (lambda.min):\n")
cat(sprintf("  R-squared: %.4f\n", r2_min))
cat(sprintf("  Intercept: %.4f (true: %.2f)\n", lasso_intercept_min, gamma_0))
cat("\n")

cat("LASSO (lambda.1se):\n")
cat(sprintf("  R-squared: %.4f\n", r2_1se))
cat(sprintf("  Intercept: %.4f (true: %.2f)\n", lasso_intercept_1se, gamma_0))
cat("\n")

cat(rep("=", 70), "\n", sep = "")

# ============================================================================
# 12. Visualizations: LASSO Coefficients and clamp_plot()
# ============================================================================

cat("\n\n")
cat(rep("=", 70), "\n", sep = "")
cat("VISUALIZATIONS\n")
cat(rep("=", 70), "\n", sep = "")
cat("\n")

# Source clamp_plot function
source("../R/clamp_plot.R")

# Plot 1: LASSO Coefficients (ranked by absolute value)
cat("Creating LASSO coefficient plots...\n")

pdf(file.path(plots_dir, "lasso_coefficients.pdf"), width = 12, height = 6)
par(mfrow = c(1, 2), mar = c(10, 4, 3, 2))

# Lambda.min
lasso_abs <- abs(lasso_coef_min)
lasso_rank <- order(lasso_abs, decreasing = TRUE)
lasso_colors <- ifelse(1:p %in% true_causal_vars, "red", "gray70")

barplot(lasso_coef_min[lasso_rank],
        names.arg = colnames(X)[lasso_rank],
        las = 2,
        col = lasso_colors[lasso_rank],
        main = "LASSO Coefficients (lambda.min)",
        ylab = "Coefficient Value",
        cex.names = 0.8)
abline(h = 0, lty = 2)
legend("topright", legend = c("True Causal", "Non-Causal"),
       fill = c("red", "gray70"), cex = 0.8)

# Lambda.1se
lasso_abs_1se <- abs(lasso_coef_1se)
lasso_rank_1se <- order(lasso_abs_1se, decreasing = TRUE)

barplot(lasso_coef_1se[lasso_rank_1se],
        names.arg = colnames(X)[lasso_rank_1se],
        las = 2,
        col = lasso_colors[lasso_rank_1se],
        main = "LASSO Coefficients (lambda.1se)",
        ylab = "Coefficient Value",
        cex.names = 0.8)
abline(h = 0, lty = 2)
legend("topright", legend = c("True Causal", "Non-Causal"),
       fill = c("red", "gray70"), cex = 0.8)

par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
dev.off()

cat(sprintf("LASSO plots saved to: %s\n\n",
            file.path(plots_dir, "lasso_coefficients.pdf")))

# Plot 2: CLAMP results using clamp_plot()
cat("Creating CLAMP plot...\n")
pdf(file.path(plots_dir, "clamp_plot.pdf"), width = 10, height = 6)
clamp_plot(fit, y = "PIP", add_legend = "topright",
           main = "CLAMP Results (with Overlap Weights)")
dev.off()
cat(sprintf("CLAMP plot saved to: %s\n\n",
            file.path(plots_dir, "clamp_plot.pdf")))

# Plot 3: SuSiE results using clamp_plot()
cat("Creating SuSiE plot using clamp_plot()...\n")
pdf(file.path(plots_dir, "susie_plot.pdf"), width = 10, height = 6)
clamp_plot(fit_susie, y = "PIP", add_legend = "topright",
           main = "SuSiE Results (Standard Method)")
dev.off()
cat(sprintf("SuSiE plot saved to: %s\n\n",
            file.path(plots_dir, "susie_plot.pdf")))

# Plot 4: Side-by-side comparison of PIPs
cat("Creating PIP comparison plot...\n")

pdf(file.path(plots_dir, "pip_comparison.pdf"), width = 12, height = 6)
par(mfrow = c(1, 1), mar = c(10, 4, 3, 2))

# Order by CLAMP PIP
pip_rank <- order(fit$pip, decreasing = TRUE)

# Create matrix for barplot
pip_matrix <- rbind(
  CLAMP = fit$pip[pip_rank],
  SuSiE = fit_susie$pip[pip_rank]
)

barplot(pip_matrix,
        beside = TRUE,
        names.arg = colnames(X)[pip_rank],
        las = 2,
        col = c("steelblue", "coral"),
        main = "Posterior Inclusion Probabilities: CLAMP vs SuSiE",
        ylab = "PIP",
        cex.names = 0.8,
        ylim = c(0, 1.1))
abline(h = 0.5, lty = 2, col = "black", lwd = 2)
legend("topright", legend = c("CLAMP", "SuSiE", "Threshold (0.5)"),
       fill = c("steelblue", "coral", NA),
       border = c("black", "black", NA),
       lty = c(NA, NA, 2),
       lwd = c(NA, NA, 2),
       cex = 0.8)

# Highlight true causal variables
true_positions <- which(pip_rank %in% true_causal_vars)
for (pos in true_positions) {
  rect((pos - 1) * 3 + 0.5, -0.05, pos * 3 + 0.5, 1.15,
       border = "red", lwd = 2)
}

par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
dev.off()

cat(sprintf("PIP comparison plot saved to: %s\n\n",
            file.path(plots_dir, "pip_comparison.pdf")))

cat(rep("=", 70), "\n", sep = "")
cat("ALL VISUALIZATIONS COMPLETE\n")
cat(rep("=", 70), "\n", sep = "")
