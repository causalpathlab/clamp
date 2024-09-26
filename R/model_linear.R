#' @title Functions related to linear regression model
#'
#' @description This file defines functions related to linear regression model,
#' which is equivalent to applying a generalized linear model with the identity
#' link function.
#'

# Expected log-likelihood for a (ordinary) susie fit for
# a linear regression model
#' @keywords internal
Eloglik_linear = function (X, Y, s) {
  n = nrow(X)
  return(-(n/2) * .loge(2*pi*s$sigma2) - (1/(2*s$sigma2)) * get_ER2(X,Y,s))
}


# Expected squared residuals.
# s$Xr is column sum of Xr_L
get_ER2 = function (X, Y, s) {
  Xr_L = compute_MXt(s$alpha * s$mu, X) # L by N matrix
  postb2 = s$alpha * s$mu2 # Posterior second moment.
  return(sum((Y - s$Xr)^2) - sum(Xr_L^2) + sum(attr(X,"d") * t(postb2)))
}
