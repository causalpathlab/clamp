#' @title Functions related to binomial model with logistic link
#'
#' @description This file defines functions related to binomial model with
#' logistic link, including (log of) pseudo-variance and pseudo-response
#' that are needed to be calculated during iterative
#' fitting of the (approximated) weighted linear regression model,
#' and the exact and approximated log-likelihood.
#'
#' The functions below output vectors of the same length of y

#' @keywords internal
#'
expit <- function(eta) {
  res <- ifelse(eta > 0,
                1 / (1 + exp(-eta)),
                exp(eta) / (1 + exp(eta)))
  return(res)
}

#' \code{compute_logw2_logistic} computes the log-pseudo-variance
#' for each observation
#' @keywords internal
#'
log_pseudo_variance_logistic <- function(eta) {
  # Input \code{y} is not applied in this function.

  logw2 <- ifelse(eta > 1e2,
                  eta,
                  2 * log1p(exp(eta)) - eta)
  return(logw2)
}

#' \code{compute_psdresponse_logistic} computes the pseudo-response
#' for each observation
#' @keywords internal
#'
pseudo_response_logistic <- function(eta, y) {

  if (length(eta) != length(y))
    stop("Dimensions of input linear predictor eta and y do not match")

  prob <- expit(eta)                    # estimated probability
  logw2 <- compute_logw2_logistic(eta)  # log of pseudo-variance
  zz <- eta + exp(logw2) * (y - prob)   # pseudo-response

  return(zz)
}

#' @keywords internal
#'
loglik_logistic <- function(eta, y) {

  if (dim(eta)[1] != length(y))
    stop("Dimensions of input linear predictor eta and y do not match")

  res <- ifelse(eta <= 500,
                y * eta - log1p(exp(eta)),
                y * eta - eta)  # for computational stability.
  return(res)
}
