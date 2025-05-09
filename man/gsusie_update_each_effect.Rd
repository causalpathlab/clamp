% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gsusie_update_each_effect.R
\name{gsusie_update_each_effect}
\alias{gsusie_update_each_effect}
\title{Update each effect once in a generalized linear model; each sub-model is a weighted SER.}
\usage{
gsusie_update_each_effect(
  X,
  y,
  s,
  model,
  estimate_prior_variance = FALSE,
  estimate_prior_method = "optim",
  check_null_threshold = 0,
  abnormal_proportion = 0.5,
  robust_method = c("none", "huber"),
  robust_estimator = c("M", "S")
)
}
\arguments{
\item{X}{An (n by p) matrix of regressor variables}

\item{y}{An n vector of response variable}

\item{s}{A gsusie fit}

\item{estimate_prior_variance}{boolean indicating whether to
estimate prior variance}

\item{check_null_threshold}{Float, a threshold on the log scale to
compare likelihood between current estimate and zero the null}

\item{abnormal_proportion}{a value between 0 and 1. If the number of detected
abnormal subjects exceeds \eqn{abnormal_proportion * nrow(X)},
stop fitting the model.}

\item{robust_method}{A string, whether and which robust method is applied
when fitting the model. \code{robust_method="none"} specifies that no
robust method is applied. \code{robust_method="huber"} specifies that the
Huber weighting method is applied.}

\item{robust_estimator}{A string, which robust estimator is applied.
\code{robust_estimator="M"} indicates the M-estimator is applied, and
\code{robust_estimator="S"} indicates the S-estimator is applied.}
}
\description{
Update each effect once in a generalized linear model; each sub-model is a weighted SER.
}
