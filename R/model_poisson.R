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
#' pseudo-response; and \code{loglik_poisson} computes the log-likelihood
#' function.
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
  clip_eta <- ifelse(eta < eta_tol, eta_tol, eta)

  return(eta + y * exp(-clip_eta) - 1)
  # return(eta + y * exp(-eta) - 1)
}


loglik_poisson <- function(eta, y, eta_tol = 50) {

  if (length(eta) != length(y))
    stop("Dimensions of input eta and y do not match")

  ## clip extremely large exp(eta)
  clip_eta <- ifelse(eta > eta_tol, eta_tol, eta)

  res <- y * eta - exp(clip_eta) + lfactorial(y)

  return(res)
}


# compute_loglik_apprx_poisson <- function(eta, y) {
#   logw2 <- compute_logw2_poisson(eta)
#   zz <- compute_psdresponse_poisson(eta, y)  # pseudo-response
#   res <- - 1 / 2 * (eta - zz)^2 * exp(-logw2)
#   return(res)
# }
