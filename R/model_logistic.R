#' @title Functions related to binomial model with logistic link
#'
#' @description This file defines functions related to binomial model with
#' logistic link, including (log of) pseudo-variance and pseudo-response
#' that are needed to be calculated during iterative
#' fitting of the (approximated) weighted linear regression model,
#' and the exact and approximated log-likelihood.
#'
#' The functions below output vectors of the same length of y
#'
#' \code{log_pseudo_variance_logistic} computes the log of pseudo-variance of
#' the logistic regression model; \code{pseudo_response_logistic} computes the
#' pseudo-response; \code{loglik_logistic} computes the log-likelihood
#' function; and \code{inverse_link_logistic} computes the inverse of link
#' function, i.e., the function of estimated expectation.
#'
#' @keywords internal
#'
expit <- function(eta) {
  res <- ifelse(eta > 0,
                1 / (1 + exp(-eta)),
                exp(eta) / (1 + exp(eta)))
  return(res)
}
#'
#' @keywords internal
#'
log_pseudo_variance_logistic <- function(eta) {
  # Input \code{y} is not applied in this function.

  logw2 <- ifelse(eta > 1e2,
                  eta,
                  2 * log1p(exp(eta)) - eta)
  return(logw2)
}

#'
#'
#' @keywords internal
#'
pseudo_response_logistic <- function(eta, y) {

  if (length(eta) != length(y))
    stop("Dimensions of input linear predictor eta and y do not match")

  prob <- expit(eta)                    # estimated probability
  logw2 <- log_pseudo_variance_logistic(eta)  # log of pseudo-variance
  zz <- eta + exp(logw2) * (y - prob)   # pseudo-response

  return(zz)
}

#' @keywords internal
#' returns an n-dim vector
loglik_logistic <- function(X, y, s) {

  ## standardize X if applicable
  X_ <- sweep(X, 2, attr(X, "scaled:center"), "-")  ## centralized
  X_ <- sweep(X_, 2, attr(X, "scaled:scale"), "/")

  # compute the linear predictor
  eta <- tcrossprod(X_, t(colSums(s$alpha * s$mu)))

  ## The following lines should be equivalent:
  ## eta <- tcrossprod(X, t(colSums(s$alpha*s$mu) / attr(X, "scaled:scale"))) -
  ##   attr(X, "scaled:center") * (colSums(s$alpha*s$mu)/attr(X, "scaled:scale"))

  # compute the log-likelihood
  res <- ifelse(eta <= 500,
                y * eta - log1p(exp(eta)),
                y * eta - eta)  # for computational stability.
  return(res)
}

#' @keywords internal
#'
inverse_link_logistic <- function(eta) {
  return(expit(eta))
}
