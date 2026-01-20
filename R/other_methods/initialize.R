#' @title Initialize a clamp object using regression coefficients
#' Copied and modified from susieR
#'
#' @param coef_index An L-vector containing the the indices of the
#'   nonzero coefficients.
#'
#' @param coef_value An L-vector containing initial coefficient
#' estimates.
#'
#' @param p A scalar giving the number of variables.
#'
#' @return A list with elements \code{alpha}, \code{mu} and \code{mu2}
#'   to be used by \code{clamp}.
#'
#' @export
#'
clamp_init_coef = function (coef_index, coef_value, p) {
  maxL = length(coef_index)
  if (L <= 0)
    stop("Need at least one non-zero effect")
  if (!all(coef_value != 0))
    stop("Input coef_value must be non-zero for all its elements")
  if (maxL != length(coef_value))
    stop("Inputs coef_index and coef_value must of the same length")
  if (max(coef_index) > p)
    stop("Input coef_index exceeds the boundary of p")
  alpha = matrix(0,nrow = maxL,ncol = p)
  mu = matrix(0,nrow = maxL,ncol = p)
  for(i in 1:L){
    alpha[i,coef_index[i]] = 1
    mu[i,coef_index[i]] = coef_value[i]
  }
  out = list(alpha = alpha,mu = mu,mu2 = mu*mu)

  class(out) = c("clamp","list")
  return(out)
}

# Set default clamp initialization.
init_setup <- function (n, p, maxL, family,
                        scaled_prior_variance, residual_variance,
                        prior_inclusion_prob, varY, standardize) {
  if (!is.numeric(scaled_prior_variance) || scaled_prior_variance < 0)
    stop("Scaled prior variance should be positive number")

  if (scaled_prior_variance > 1 && standardize)
    stop("Scaled prior variance should be no greater than 1 when ",
         "standardize = TRUE")

  if(is.null(residual_variance)) {

    if (family == "linear") {
      ## for linear models,
      ## 1. initialize residual_variance as varY
      if (is.null(residual_variance)) {
        residual_variance = varY
      }
      ## 2. initialize prior_varD as follows...
      prior_varD = scaled_prior_variance * varY
    } else {
      ## for logistic and poisson regression models
      ## for generalized linear models,
      ## 1. initialize residual_variance as 1
      residual_variance = 1
      ## 2. initialize prior_varD as 1
      prior_varD = 1
    }

  }

  if(is.null(prior_inclusion_prob)){
    prior_inclusion_prob = rep(1/p,p)
  } else {
    if(all(prior_inclusion_prob == 0)){
      stop("Prior weight should greater than 0 for at least one variable.")
    }
    prior_inclusion_prob = prior_inclusion_prob / sum(prior_inclusion_prob)
  }
  if(length(prior_inclusion_prob) != p)
    stop("Prior weights must have length p")
  if (p < maxL) maxL = p

  s = list(family = family,
           alpha  = matrix(1/p,nrow = maxL,ncol = p),
           mu     = matrix(0,nrow = maxL,ncol = p),
           mu2    = matrix(0,nrow = maxL,ncol = p),
           deltahat = matrix(0, nrow = maxL, ncol = p),  ## MLE
           shat2 = matrix(0, nrow = maxL, ncol = p), ## shat2 of MLE
           Xr     = rep(0,n),
           # KL     = rep(as.numeric(NA),maxL),
           # EWR2 = rep(as.numeric(NA), maxL),
           EWR2 = rep(as.numeric(NA), maxL),
           Eloglik = rep(as.numeric(NA), maxL),
           logBF  = rep(as.numeric(NA),maxL),
           logBF_variable = matrix(as.numeric(NA),maxL,p),
           sigma2 = residual_variance,           ## residual variance
           prior_varD = prior_varD,         ## prior variance of coefficients b
           pie    = prior_inclusion_prob,
           abnormal_subjects = NULL,        ## abnormal subject indices
           importance_weight = rep(1, n))   ## importance weights for robust methods

  if (family == "linear") class(s) = "clamp"
  else class(s) = "gsusie"

  return(s)
}

# Update a clamp fit object in order to initialize clamp model.
init_finalize = function (s, X = NULL, Xr = NULL) {
  if(length(s$prior_varD) == 1)
    s$prior_varD = rep(s$prior_varD, nrow(s$alpha))

  # Check sigma2.
  if (!is.numeric(s$sigma2))
    stop("Input residual variance sigma2 must be numeric")

  # Avoid problems with dimension if input is a 1 x 1 matrix.
  s$sigma2 = as.numeric(s$sigma2)
  if (length(s$sigma2) != 1)
    stop("Input residual variance sigma2 must be a scalar")
  if (s$sigma2 <= 0)
    stop("Residual variance sigma2 must be positive (is your var(Y) zero?)")

  # check prior variance
  if (!is.numeric(s$prior_varD))
    stop("Input prior variance must be numeric")
  if (!all(s$prior_varD >= 0))
    stop("prior variance must be non-negative")
  if (!all(dim(s$mu) == dim(s$mu2)))
    stop("dimension of mu and mu2 in input object do not match")
  if (!all(dim(s$mu) == dim(s$alpha)))
    stop("dimension of mu and alpha in input object do not match")
  if (nrow(s$alpha) != length(s$prior_varD))
    stop("Input prior variance prior_varD must have length of nrow of alpha in ",
         "input object")

  # Update Xr.
  if (!missing(Xr))
    s$Xr = Xr
  if (!missing(X))
    s$Xr = compute_Xb(X,colSums(s$mu * s$alpha))

  # Reset KL and lbf.
  # s$KL = rep(as.numeric(NA),nrow(s$alpha))
  s$logBF = rep(as.numeric(NA),nrow(s$alpha))
  class(s) = "clamp"
  return(s)
}
