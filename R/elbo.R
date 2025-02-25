# @title Get objective function from data and clamp fit object.
# @param data A flash data object.
# @param f A flash fit object.
# @keywords internal

# get_elbo = function (X, y, s, W=NULL, model) {
get_elbo = function(s) {

  post_varB <- s$mu2 - s$mu^2  # avoid 0

  revKL <- reverseKL_va_vs_prior(alpha = s$alpha,
                                 post_varB = post_varB,
                                 mu = s$mu,
                                 pie = s$pie,
                                 prior_varB = s$prior_varB)

  Eloglik <- s$Eloglik

  # return(sum(model$loglik(X,Y,s)) + sum(revKL))
  return( mean(Eloglik) + sum(revKL) )
}

# to avoid log(0)...
.loge <- function(t) {
  log(t+.Machine$double.eps)
}

## Reverse KL of the variational and prior distributions of all layers
#' @param alpha     an L by p matrix, posterior inclusion probabilities
#' @param post_varB an L by p matrix, posterior variance of Beta
#' @param mu        an L by p matrix, posterior mean of Beta
#' @param pie       a  p-dim vector, prior inclusion probabilities
#' @param prior_varB an L-dim vector, prior variance of Beta (one for each layer)
#'
reverseKL_va_vs_prior <- function(alpha, post_varB, mu, pie, prior_varB) {

  rowSums(
    0.5*alpha*( 1 + sweep(.loge(post_varB), 1, .loge(prior_varB), "-") -
                  sweep(mu^2+post_varB, 1, prior_varB+.Machine$double.eps, "/")  ) +
      alpha * sweep(-.loge(alpha), 2, .loge(pie), "+")
  )
}

## If variances of variational (posterior) and prior distribution are both 0,
## then this function will return a certain fixed value (0.5).

