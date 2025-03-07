#' @rdname eloglik_wser
#' @title Functions related to IPW single effect regression model
#'
#' @description This file defines functions related to weighted linear
#' regression model
#'

# Expected (composite) log-likelihood of layer-l weighted single effect regression.
#'
#' @param X (n by p) matrix
#' @param y (n by 1) vector
#' @param W (n by p) matrix of weights
#' @param s a clamp fit, containing a list of posterior distributions
#' @param alpha a p-vector of posterior inclusion probabilities
#' @param mu a p-vector of posterior means (first moments)
#' @param mu2 a p-vector of posterior second moments
#'
#' @keywords internal
WSER_posterior_e_loglik = function (X, y, W, alpha, mu, mu2,
                         residual_variance) {
  n = nrow(X)

  EWR2 <- get_EWR2_l(X, y, W, alpha, mu, mu2)

  return(-n/2*log(2*pi*residual_variance) + 1/2*sum(alpha*colSums(log(W))) -
           1/(2*residual_variance) * EWR2)
}

#' Expectation of weighted sum of squared residual (WRSS) of layer l
#'
#'  X: (n by p) matrix
#'  y: (n by 1) vector
#'  W: (n by p) matrix of weights
#'  alpha: a p-vector of posterior inclusion probabilities
#'  mu: a p-vector of posterior means (first moments)
#'  mu2: a p-vector of posterior second moments
#'
get_EWR2_l <- function(X, y, W,
                       alpha=NULL, mu=NULL, mu2=NULL) {

  if (is.null(alpha) | is.null(mu) | is.null(mu2)) {
    stop("Specify `alpha, mu, mu2`!")
  }

  X_ <- sapply(1:ncol(X), function(j) {X[,j] - weighted.mean(x=X[,j], w=W[,j])})
  y_ <- sapply(1:ncol(X), function(j) {y - weighted.mean(x=y, w=W[,j])})

  postb  <- alpha * mu
  postb2 <- alpha * mu2

  EWR2 <- sum(sweep(W*(y_^2), 2, alpha, "*")) -
    2 * sum(sweep(W*X_*y_, 2, postb, "*")) +
    sum(sweep(W*(X_^2), 2, postb2, "*"))

  return(EWR2)
}
