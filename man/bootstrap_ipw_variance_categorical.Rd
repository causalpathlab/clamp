% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootstrap_ipw_variance_categorical.R
\name{bootstrap_ipw_variance_categorical}
\alias{bootstrap_ipw_variance_categorical}
\title{Bootstrap the variance of inverse probability weighting (IPW) estimator}
\usage{
bootstrap_ipw_variance_categorical(
  X,
  y,
  W,
  causal_effect_estimator = c("mHT"),
  nboots = 100,
  seed = NULL,
  baseline_bootMeans = TRUE
)
}
\arguments{
\item{X}{An n by p data matrix}

\item{y}{An n-dim vector}

\item{W}{An n by p weight matrix, each entry \eqn{w_{ij}} is the inverse of
propensity score.}

\item{nboots}{The number of bootstrap replicates (samples). By default,
\code{nboots=100}.}

\item{seed}{Random seed. If \code{seed=NULL}, its default is \code{Sys.time()}.}

\item{baseline_bootMeans}{Whether to compute the bootstrap means of the
baseline levels}
}
\value{
The function returns a p-dim vector of estimated variances,
each entry of which is the bootstrap variance of \eqn{X_{.j}}.
}
\description{
This function estimates the variance of each regression coefficient,
by (nonparametric) bootstrapping.
}
\keyword{internal}
