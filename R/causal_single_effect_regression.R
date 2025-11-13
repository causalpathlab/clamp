#' @rdname causal_single_effect_regresion
#'
#' @title Single effect regression with inverse probability weighting estimates for categorical treatments
#'
#' @description
#' This function is to compute the posterior distribution of the regression
#' coefficients of a Causal-SER model. Adapted from:
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
#' @param causal_effect_estimator The estimation method of the causal effect.
#' \code{"ipw"} estimates the ATE with inverse probability weighting estimator.
#'
#' @param variance_estimator The estimation method of the variance of the
#' (Horvitz-Thompson) ATE estimators.
#' \code{"bootstrap"} estimates the variances from \code{nboots} bootstrap replicates.
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
#' \item{deltahat}{Vector of (non Bayesian) average treatment estimates (ATEs).}
#'
#' \item{shat2}{Vector of approximated variance estimates of ATEs}
#'
#' \item{thetahat0}{Vector of average causal effect estimates of baselines.}
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
causal_single_effect_regression <-
  function(y, X,
           W = NULL,
           prior_varD,
           residual_variance = 1,
           prior_inclusion_prob = NULL,
           optimize_prior_varD = c("none", "optim", "EM", "simple"),
           check_null_threshold = 0,
           causal_effect_estimator = c("mHT"),
           variance_estimator = c("bootstrap"),
           nboots = 100,
           seed = NULL) {

    optimize_prior_varD <- match.arg(optimize_prior_varD)

    nn <- nrow(X)

    # Check W
    if (is.null(W)) {
      W <- rep(1, times = nn)
    } else if (!all(dim(W) == dim(X))) {  ## (n by p) matrix
      stop("Dimensions of W and X do not match!")
    }

    switch(variance_estimator,
           "bootstrap" = {
             ATE <-
               estimate_average_treatment_effect_categorical(X, y, W,
                    causal_effect_estimator = causal_effect_estimator)
             deltahat <- ATE$deltahat
             # thetahat0 <- ATE$thetahat0

             boot_result <-
               bootstrap_ipw_variance_categorical(X = X, y = y, W = W,
                            causal_effect_estimator = causal_effect_estimator,
                            nboots = nboots, seed = seed)
             shat2 <- boot_result$deltahat_bootVars
             # thetahat0_bootMeans <- boot_result$thetahat0_bootMeans
           }
           )

    plevels <- length(deltahat)  ## number of levels
    # plevels <- length(attr(X, "K_minus_1_dummy_indices"))

    # Check prior_inclusion_probability
    ## If `prior_inclusion_probability` is consistent to all variables,
    ## its value should not affect the posterior inclusion probability (alpha).
    if (is.null(prior_inclusion_prob)) {
      prior_inclusion_prob <- rep(1/plevels, times = plevels)
    } else {
      if (length(prior_inclusion_prob) != plevels)
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
    posterior_varD <- 1 / (1/shat2 + 1/prior_varD)   # posterior variance
    mu <- prior_varD / (prior_varD + shat2) * deltahat  # posterior mean
    mu2 <- posterior_varD + mu^2                      # posterior second moment

    # ABF for WSER model
    logBF_model <- maxlogPO + log(sum(logPO_weighted))
    # = log(sum(ABF x prior_weights))

    if (optimize_prior_varD == "EM") {
      prior_varD <- optimize_prior_variance(optimize_prior_varD, deltahat, shat2,
                                            prior_inclusion_prob,
                                            alpha, mu2,
                                    check_null_threshold = check_null_threshold)
    }

    return(list(alpha = alpha,      # level-wise posterior inclusion probability
                mu = mu,            # posterior mean
                mu2 = mu2,          # posterior second moment
                deltahat = deltahat,  # (non bayesian) ATEs
                shat2 = shat2,        # approximate variances of ATEs
                # thetahat0 = thetahat0,  # (non bayesian) baseline causal effect estimates. Replaced by the bootstrap means?
                logBF = logBF,      # layer-wise log of Bayes factor
                logBF_model = logBF_model,  # log of Bayes factor of model
                prior_varD = prior_varD  # prior variance of coefficients B
    ))
  }


