# @title Get objective function from data and clamp fit object.
# @param data A flash data object.
# @param f A flash fit object.
# @keywords internal
get_elbo = function (X, Y, s, W=NULL, model) {
  include_idx <- s$prior_varB > 1e-9

  post_varB <- s$mu2 - s$mu^2  # avoid 0
  revKL <- sapply(1:nrow(s$alpha),
                  function(l){
                    reverseKL_l_fun(alpha.l      = s$alpha[l,],
                                    post_varB.l  = post_varB[l,],
                                    mu.l         = s$mu[l,],
                                    pie.l        = s$pie[l],
                                    prior_varB.l = s$prior_varB[l])})
  revKL <- revKL[include_idx]

  return(sum(model$loglik(X,Y,s)) + sum(revKL))  ##Eloglik!!
}

# to avoid log(0)...
.loge <- function(t) {
  log(t+.Machine$double.eps)
}

# Compute the (layer-wise) reverse KL of the posterior and prior distribution
# @param alpha.l     an p-dim vector, posterior inclusion probabilities
# @param post_varB.l an p-dim vector, posterior variance of B
# @param mu.l         a p-dim vector, poterior mean of B
# @param pie.l        a number, prior inclusion probabilities
# @param prior_varB.l a number, prior variance of B
reverseKL_l_fun <- function(alpha.l, post_varB.l, mu.l,
                                 pie.l, prior_varB.l) {
  sum(0.5*alpha.l * (1+.loge(post_varB.l)-.loge(prior_varB.l) -
                         (mu.l^2+post_varB.l) / prior_varB.l) +
        alpha.l * (.loge(pie.l) - .loge(alpha.l)))
}

# @title posterior expected log-likelihood for a single effect regression
# @param X an n by p matrix of covariates
# @param Y an n vector of regression outcome
# @param s2 the residual variance
# @param Eb the posterior mean of b (p vector) (alpha * mu)
# @param Eb2 the posterior second moment of b (p vector) (alpha * mu2)
SER_posterior_e_loglik = function (X, Y, s2, Eb, Eb2) {
  n = nrow(X)
  return(-0.5*n*log(2*pi*s2) - 0.5/s2*(sum(Y*Y)
                                       - 2*sum(Y*compute_Xb(X,Eb))
                                       + sum(attr(X,"d") * Eb2)))
}



