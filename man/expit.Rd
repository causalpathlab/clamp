% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_logistic.R
\name{expit}
\alias{expit}
\title{Functions related to binomial model with logistic link}
\usage{
expit(eta)
}
\description{
This file defines functions related to binomial model with
logistic link, including (log of) pseudo-variance and pseudo-response
that are needed to be calculated during iterative
fitting of the (approximated) weighted linear regression model,
and the exact and approximated log-likelihood.

The functions below output vectors of the same length of y

\code{log_pseudo_variance_logistic} computes the log of pseudo-variance of
the logistic regression model; \code{pseudo_response_logistic} computes the
pseudo-response; \code{loglik_logistic} computes the log-likelihood
function; and \code{inverse_link_logistic} computes the inverse of link
function, i.e., the function of estimated expectation.
}
\keyword{internal}
