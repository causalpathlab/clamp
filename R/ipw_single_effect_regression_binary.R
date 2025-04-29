#' @rdname ipw_single_effect_regresion_binary
#'
#' @title Single effect regression with inverse probability weighting estimates for binary treatments
#'
#' @description
#' This function is to compute the posterior distribution of the regression
#' coefficients of a IPW-SER model. Adapted from:
#' <https://github.com/stephenslab/susieR/blob/master/R/single_effect_regression.R>
#' This is basically for function \code{clamp}.
#'
#' @param y An (n by 1) vector of responses.
#'
#' @param X An (n by p) data matrix; each entry would be 0 or 1.
#'
#' @param W An (n by p) matrix; each entry is the weight of each each
#' entry of \eqn{X}; it should be of the same size as \eqn{X}.
#'
#' @param residual_variance A scalar refers to the residual variance.
#' To be consistent with susieR, it should be specified when \code{W=NULL}.
#'
#' @param prior_inclusions_prob An (p by 1) vector of prior inclusion
#' probabilities.
#'
#' @param prior_varD A scalar given the (initial) prior variance of the coefficient.
#'
#' @param optimize_prior_varD The optimization method to use for fitting the prior
#' variance.
#'
#' @param check_null_threshold Scalar specifying the threshold on the log-scale
#' to compare likelihood between current estimate and zero (null).
#'
#' @param mle_estimator The estimation method of the MLEs of the ATEs.
#' \code{"ipw"} estimates the ATE with inverse probability weighting estimator,
#' and \code{"wls"} applies the equivalent weighted least-squares estimator.
#' The estimates from these two values should be the same if the input X is not
#' centralized and standardized.
#'
#' @param mle_variance_estimator The estimation method of the variance of the MLEs,
#' or equivalently, the variance of (Horvitz-Thompson) ATE estimators.
#' \code{"bootstrap"} estimates the variances from \code{nboots} bootstrap replicates,
#' \code{"sandwich"} uses sandwich robust variance estimators.
#' \code{"naive"} calculates the variances from the weighted least square model,
#' treating weights as constants (not recommended).
#'
#' @param nboots The number of bootstrap replicates. By default, \code{nboots=100}.
#'
#' @param seed Random seed utilized when \code{mle_variance_estimator="bootstrap"}.
#' If \code{seed=NULL}, its default is \code{Sys.time()}.
#'
#' @returns A list with the following elements:
#'
#' \item{alpha}{Vector of posterior inclusion probabilities;
#'  \code{alpha[i]} is posterior probability that the ith
#'  coefficient is non-zero.}
#'
#' \item{mu}{Vector of posterior means (conditional on inclusion).}
#'
#' \item{mu2}{Vector of posterior second moments (conditional on inclusion).}
#'
#' \item{logBF}{Vector of log of (asymptotic) Bayes factor for each variable.}
#'
#' \item{logBF_model}{Log of (asymptotic) Bayes factor for the
#'          weighted single effect regression.}
#'
#' \item{prior_varD}{Prior variance (after optimization if \code{optimize_prior_varD != "none"}).}
#'
#' @importFrom stats dnorm
#' @importFrom stats uniroot
#' @importFrom stats optim
#' @importFrom Matrix colSums
#' @importFrom sandwich vcovHC
#'
#' @keywords internal
#'
ipw_single_effect_regression_binary <-
  function(y, X,
           W = NULL,
           prior_varD,
           residual_variance = 1,
           prior_inclusion_prob = NULL,
           optimize_prior_varD = c("none", "optim", "EM", "simple"),
           check_null_threshold = 0,
           mle_estimator = c("mHT", "WLS"),
           mle_variance_estimator = c("bootstrap", "sandwich"),
           nboots = 100,
           seed = NULL) {

    optimize_prior_varD <- match.arg(optimize_prior_varD)

    n <- nrow(X)
    p <- ncol(X)

    # Check W
    if (is.null(W)) {
      W <- rep(1, times = n)
    } else if (!all(dim(W) == dim(X))) {  ## (n by p) matrix
      stop("Dimensions of W and X do not match!")
    }

    switch(mle_variance_estimator,
           "bootstrap" = {
             deltahat <- estimate_average_treatment_effect_binary(X, y, W,
                                                  mle_estimator = mle_estimator)
             shat2 <-
               bootstrap_ipw_variance_binary(X = X, y = y, W = W,
                                      mle_estimator = mle_estimator,
                                      nboots = nboots, seed = seed)
           },

           "sandwich" = {
             # Both X and y need to be centralized!
             deltahat <- rep(NA, times = p)
             shat2 <- rep(NA, times = p)
             for (j in 1 : p) {
               lmfit <- lm(y ~ X[,j], weights = W[,j])
               deltahat[j] <- coefficients(lmfit)[2]
               shat2[j] <- vcovHC(lmfit, type = "HC")[2,2]
             }
           })

    # Check prior_inclusion_probability
    ## If `prior_inclusion_probability` is consistent to all variables,
    ## its value should not affect the posterior inclusion probability (alpha).
    if (is.null(prior_inclusion_prob)) {
      prior_inclusion_prob <- rep(1/p, times = p)
    } else {
      if (length(prior_inclusion_prob) != p)
        stop("The length of prior inclusion probability does not match")
    }

    if (optimize_prior_varD != "EM" && optimize_prior_varD != "none") {
      prior_varD <- optimize_prior_variance(optimize_prior_varD, deltahat, shat2,
                                  prior_inclusion_prob,
                                  alpha = NULL, mu2 = NULL,
                                  prior_varD_init = prior_varD,
                                  check_null_threshold = check_null_threshold)
    }

    # compute log of Bayes factor (logBF) and log of posterior odds (logPO)
    # .loge <- function(t) log(t+.Machine$double.eps) (in `elbo.R`)
    zscore2 <- deltahat^2 / shat2
    logBF <- 1/2 * (.loge(shat2 /(shat2 + prior_varD)) +
                      zscore2 * prior_varD / (prior_varD + shat2))
    logPO <- logBF + log(prior_inclusion_prob + sqrt(.Machine$double.eps))
    # deal with special case of infinite shat2
    # (e.g. happens if X does not vary)
    logBF[is.infinite(shat2)] <- 0
    logPO[is.infinite(shat2)] <- 0
    maxlogPO <- max(logPO)

    # logPO_tilde is proportional to posterior odds = BF * prior,
    # but subtract max for numerical stability
    logPO_weighted <- exp(logPO - maxlogPO)

    # Update the posterior estimates
    # Posterior prob for each variable
    alpha <- logPO_weighted / sum(logPO_weighted)    # posterior inclusion probability
    post_varD <- 1 / (1/shat2 + 1/prior_varD)        # posterior variance
    mu <- prior_varD / (prior_varD + shat2) * deltahat  # posterior mean
    mu2 <- post_varD + mu^2                      # posterior second moment

    # ABF for WSER model
    logBF_model <- maxlogPO + log(sum(logPO_weighted))
    # = log(sum(ABF x prior_weights))

    if (optimize_prior_varD == "EM") {
      prior_varD <- optimize_prior_variance(optimize_prior_varD, deltahat, shat2,
                                            prior_inclusion_prob,
                                            alpha, mu2,
                                    check_null_threshold = check_null_threshold)
    }

    # Expected weighted sum of squared residual (EWRSS) ##??!!
    EWR2 <- get_EWR2_l(X=X, y=y, W=W,
                       alpha=alpha, mu=mu, mu2=mu2)
    # Expected log-likelihood under variational distribution ql
    Eloglik <- WSER_posterior_e_loglik(X=X, y=y, W=W,
                                       alpha=alpha, mu=mu,
                                       mu2=mu2,
                                       residual_variance=residual_variance)

    return(list(alpha = alpha,      # posterior inclusion probability
                mu = mu,            # posterior mean
                mu2 = mu2,          # posterior second moment
                deltahat = deltahat,  # maximum likelihood estimator (MLE)
                shat2 = shat2,      # approximate variance of MLE
                logBF = logBF,      # layer-wise log of Bayes factor
                logBF_model = logBF_model,  # log of Bayes factor of model
                prior_varD = prior_varD,  # prior variance of coefficients B
                EWR2 = EWR2,        # expected weighted sum of squared residual
                Eloglik = Eloglik     # expected log-likelihood
    ))
  }


