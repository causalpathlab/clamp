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
#'
#' @returns The function returns a p-dim vector of estimated variances,
#' each entry of which is the bootstrap variance of \eqn{X_{.j}}.
#'
#' @importFrom matrixStats colVars
#' @keywords internal
bootstrap_ipw_variance <- function(X, y, W, mle_estimator = c("mHT", "WLS"),
                                   nboots = 100, seed = NULL) {

  if (is.null(seed)) seed <- Sys.time()
  set.seed(seed)

  n <- nrow(X)
  p <- ncol(X)

  # Check the dimensions of W
  if (!( (length(W) == nrow(X)) | all(dim(W) == dim(X)) ) )
    stop("Dimensions of X and W do not match.")
  W <- as.matrix(W)

  # Check if each column of X is standardized or not.
  if ( all(attr(X, "scaled:scale") == 1) ) standardize=F
  else standardize = T

  # Check if each column of X is centralized or not.
  if ( all(attr(X, "scaled:center") == 0) ) centralize = F
  else centralize = T

  # bootstrap estimators
  boot_betahat <- matrix(NA_real_, nrow = nboots, ncol = p)

  for (B in 1 : nboots) {

    ind <- sample.int(n, size = n, replace = T)
    Xboot <- X[ind, , drop = F]
    yboot <- y[ind]
    Wboot <- W[ind, , drop = F]

    boot_betahat[B,] <-
      estimate_average_treatment_effect(X=Xboot, y=yboot, W=Wboot,
                                        mle_estimator=mle_estimator,
                                        centralize = centralize,
                                        standardize = standardize)
  }
  # Convert back to the original scale (if X is standardized)
  boot_betahat <- sweep(boot_betahat, 2, attr(X, "scaled:scale"), "/")

  return( colVars(boot_betahat) )
}
