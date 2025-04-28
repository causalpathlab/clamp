# @title Get objective function from data and clamp fit object.
# @keywords internal

#' @importFrom matrixStats colSums2
#' @importFrom Matrix tcrossprod
#' @importFrom Matrix crossprod

# @param X an n by p-levels matrix of covariates (K-1 dummy encoding version)
get_elbo = function(X, y, s) {

  revKL <- reverseKL_va_vs_prior(alpha = s$alpha,
                                 post_varD = s$mu2 - s$mu^2,
                                 mu = s$mu,
                                 pie = s$pie,
                                 prior_varD = s$prior_varD)

  return( Eloglik(X, y, s) + sum(revKL) )
}

# to avoid log(0)...
.loge <- function(t) {
  log(t+.Machine$double.eps)
}

## Reverse KL of the variational and prior distributions of all layers
#' @param alpha     an L by p matrix, posterior inclusion probabilities
#' @param post_varD an L by p matrix, posterior variance of Beta
#' @param mu        an L by p matrix, posterior mean of Beta
#' @param pie       a  p-dim vector, prior inclusion probabilities
#' @param prior_varD an L-dim vector, prior variance of Beta (one for each layer)
#'
reverseKL_va_vs_prior <- function(alpha, post_varD, mu, pie, prior_varD) {

  rowSums(
    0.5*alpha*( 1 + sweep(.loge(post_varD), 1, .loge(prior_varD), "-") -
                sweep(mu^2+post_varD, 1, prior_varD+.Machine$double.eps, "/")) +
      alpha * sweep(-.loge(alpha), 2, .loge(pie), "+")
  )
}
## If variances of variational (posterior) and prior distribution are both 0,
## then this function will return a certain fixed value (0.5).

# Expected loglikelihood
Eloglik <- function(X, y, s) {
  n <- nrow(X)
  return(- (n/2) * .loge(2*pi*s$sigma2) - 1/(2*s$sigma2) * get_ERSS(X, y, s))
}

# Expected residual sum of squares
# input X should be the K minus 1 dummy encoding version
get_ERSS <- function(X, y, s) {

  if (is.null(s$variable_alpha)) { ## not hieraracical PIP approach.
    postd  <- s$alpha * s$mu   ## n by pcoef matrix
    postd2 <- s$alpha * s$mu2  ## n by pcoef matrix
  } else {
    postd  <- s$variable_alpha * s$mu
    postd2 <- s$variable_alpha * s$mu2
  }

  K_minus_1_dummy_indices <- attr(X, "K_minus_1_dummy_indices")
  if (ncol(postd) != length(K_minus_1_dummy_indices)) {
    stop("Column dimension of X does not match!")
  }

  res <- sum((y - s$Xr - s$intercept)^2) +
    sum( tcrossprod(colSums2(X[, K_minus_1_dummy_indices, drop=F]^2), postd2) )-
    sum( (tcrossprod(X[, K_minus_1_dummy_indices, drop=F], postd))^2 )

  return(res)
}
