#' @title Bootstrap the variance of inverse probability weighting (IPW) estimator
#' @description This function estimates the variance of each regression coefficient,
#' by (nonparametric) bootstrapping.
#'
#' @param X An n by p data matrix
#' @param y An n-dim vector
#' @param W An n by p weight matrix, each entry \eqn{w_{ij}} is the inverse of
#' propensity score.
#' @param nboots The number of bootstrap replicates (samples). By default,
#' \code{nboots=100}.
#' @param seed Random seed. If \code{seed=NULL}, its default is \code{Sys.time()}.
#' @param baseline_bootMeans Whether to compute the bootstrap means of the
#' baseline levels
#'
#' @returns The function returns a p-dim vector of estimated variances,
#' each entry of which is the bootstrap variance of \eqn{X_{.j}}.
#'
#' @importFrom matrixStats colVars
#' @keywords internal
bootstrap_ipw_variance_categorical <- function(X, y, W,
                                   causal_effect_estimator = c("mHT"),
                                   nboots = 100, seed = NULL,
                                   baseline_bootMeans = TRUE ) {

  if (is.null(seed)) seed <- Sys.time()
  set.seed(seed)

  n <- nrow(X)

  # Check the dimensions of W
  if (!( (length(W) == nrow(X)) | all(dim(W) == dim(X)) ) )
    stop("Dimensions of X and W do not match.")
  W <- as.matrix(W)

  # bootstrap treatment effect estimates
  # ncol = number of non-reference levels in total
  boot_deltahat <- matrix(NA, nrow = nboots,
                          ncol = sum(!grepl("_0", colnames(X))))

  if (baseline_bootMeans) {
    # bootstrap baseline causal effects
    boot_thetahat0 <- matrix(NA, nrow = nboots,
                             ncol = sum(grepl("_0", colnames(X))))
  }

  for (B in 1 : nboots) {
    ind <- sample.int(n, size = n, replace = T)

    while ( any(colSums(X[ind, , drop=F]) <= 1)) {
      # if a level contains no more than 1 sample, then resample.
      # print(paste(B, ": a layer contains too little samples"))
      ind <- sample.int(n, size = n, replace = T)
    }

    Xboot <- X[ind, , drop = F]
    yboot <- y[ind]
    Wboot <- W[ind, , drop = F]

    est <- estimate_average_treatment_effect_categorical(X=Xboot, y=yboot,
                                                         W=Wboot,
                            causal_effect_estimator=causal_effect_estimator)
    boot_deltahat[B,] <- est$deltahat

    if (baseline_bootMeans) boot_thetahat0[B,] <- est$thetahat0

  }

  if ( any(is.nan(boot_deltahat)) ) {
    warning(
      sprintf("Columns include NAN bootstrap estimates: %i. \
      Variances are estimated after dropping NANs.",
              which( is.na(colVars(boot_deltahat))) ) )
  }

  colnames(boot_deltahat) <- names(est$deltahat)
  if (!baseline_bootMeans) {
    return ( colVars(boot_deltahat, na.rm=TRUE) )

  } else {
    colnames(boot_thetahat0) <- names(est$thetahat0)
    return(
      list(
        deltahat_bootVars = colVars(boot_deltahat, na.rm = TRUE),
        thetahat0_bootMeans = colMeans(boot_thetahat0, na.rm = TRUE)
        )
      )
  }

}

# Convert back to the original scale (if X is standardized)
# boot_deltahat <- sweep(boot_deltahat, 2, attr(X, "scaled:scale"), "/")

# # Check if each column of X is standardized or not.
# if ( all(attr(X, "scaled:scale") == 1) ) standardize=F
# else standardize = T
#
# # Check if each column of X is centralized or not.
# if ( all(attr(X, "scaled:center") == 0) ) centralize = F
# else centralize = T

# boot_thetahat[B,] <-
#   estimate_average_causal_effect(X=Xboot, y=yboot, W=Wboot,
#                             causal_effect_estimator=causal_effect_estimator,
#                             compute_treatment_effect = FALSE,
#                             centralize = centralize,
#                             standardize = standardize)
