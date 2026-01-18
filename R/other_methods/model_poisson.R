#' @title Functions related to Poisson regression model with log link
#'
#' @description This file defines functions related to Poisson regression model
#' with logistic link, including (log of) pseudo-variance and pseudo-response
#' that are needed to be calculated during iterative fitting of the
#' (approximated) weighted linear regression model, and the exact and
#' approximated log-likelihood.
#'
#' The functions below output vectors of the same length of y
#'
#' \code{log_pseudo_variance_poisson} computes the log of pseudo-variance of
#' the Poisson regression model; \code{pseudo_response_poisson} computes the
#' pseudo-response; \code{loglik_poisson} computes the log-likelihood
#' function; and \code{inverse_link_poisson} computes the inverse of link
#' function, i.e., the function of estimated expectation.
#'
#'
#' @keywords internal
#'
log_pseudo_variance_poisson <- function(eta) {
  return(-eta)
}

pseudo_response_poisson <- function(eta, y, eta_tol = -50) {
  if (length(eta) != length(y))
    stop("Dimensions of input eta and y do not match")

  ## clip extremely large exp(eta)
  clipped_eta <- ifelse(eta < eta_tol, eta_tol, eta)

  return(eta + y * exp(-clipped_eta) - 1)
  # return(eta + y * exp(-eta) - 1)
}


#' returns an n-dim vector
loglik_poisson <- function(X, y, s, eta_tol = 50){

  # standardize X if applicable
  X_ <- sweep(X,  2, attr(X, "scaled:center"), "-")  ## centralized
  X_ <- sweep(X_, 2, attr(X, "scaled:scale"), "/")
  # compute the linear predictor
  eta <- tcrossprod(X_, t(colSums(s$alpha * s$mu)))

  # The following lines should be equivalent:
  # eta <- tcrossprod(X, t(colSums(s$alpha*s$mu) / attr(X, "scaled:scale"))) -
  #   attr(X, "scaled:center") * (colSums(s$alpha*s$mu)/attr(X, "scaled:scale"))


  ## clip extremely large exp(eta)
  clipped_eta <- ifelse(eta > eta_tol, eta_tol, eta)
  # compute the log-likelihood
  res <- y * eta - exp(clipped_eta) - lfactorial(y)

  return(res)
}

inverse_link_poisson <- function(eta, eta_tol=50) {
  # ## clip extremely large exp(eta)
  # clipped_eta <- ifelse(eta > eta_tol, eta_tol, eta)

  return(exp(eta))
}


