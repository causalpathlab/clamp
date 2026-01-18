#' @title Functions related to linear regression model
#'
#' @description This file defines functions related to linear regression model,
#' which is equivalent to applying a generalized linear model with the identity
#' link function.
#'

# Expected log-likelihood for a (ordinary) susie fit for
# a linear regression model as a whole.
#' @keywords internal
Eloglik_linear = function (X, Y, s) {
  n = nrow(X)
  return(-(n/2) * .loge(2*pi*s$sigma2) - (1/(2*s$sigma2)) * get_ER2(X,Y,s))
}


# Expected squared residuals.
# s$Xr is column sum of Xr_L
# For family %in% c("logistic", "poisson"), there should not be a sigma2.

get_ER2 = function (X, Y, s) {

  Xr_L = compute_MXt(s$alpha * s$mu, X) # L by N matrix
  Eb2 = s$alpha * s$mu2 # Posterior second moment.

  Xd = colSums(X^2)

  # return(sum((Y - s$Xr)^2) - sum(Xr_L^2) + sum(attr(X,"d") * t(Eb2)))
  return(sum((Y - s$Xr)^2) - sum(Xr_L^2) + sum(Xd * t(Eb2)))
}


# # @title posterior expected log-likelihood for a single effect regression
# # @param X an n by p matrix of covariates
# # @param Y an n vector of regression outcome
# # @param s2 the residual variance
# # @param Eb the posterior mean of b (p vector) (alpha * mu)
# # @param Eb2 the posterior second moment of b (p vector) (alpha * mu2)
# SER_posterior_e_loglik = function (X, Y, s2, Eb, Eb2) {
#   n = nrow(X)
#   return(-0.5*n*log(2*pi*s2) - 0.5/s2*(sum(Y*Y)
#                                        - 2*sum(Y*compute_Xb(X,Eb))
#                                        + sum(attr(X,"d") * Eb2)))
# }
