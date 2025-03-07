#' @rdname weighted_single_effect_regresion
#'
#' @title Weighted single effect regression
#'
#' @description
#' The WSER function is to compute the posterior distribution of the regression
#' coefficients of a WSER model. Reference:
#' <https://github.com/stephenslab/susieR/blob/master/R/single_effect_regression.R>
#' This is basically used for function \code{gsusie}.
#'
#' @param y an (n by 1) vector of responses.
#'
#' @param X An (n by p) data matrix.
#'
#' @param W An (n by 1) vector.
#' If \code{W} is an (n by 1) vector, it contains the W of each subject;
#' for a generalized linear model using iterative reweighted least squared
#' approach, \eqn{W = exp(-logw2)};
#'
#' @param residual_variance A scalar refers to the residual variance.
#' To be consistent with susieR, it should be specified when \code{W=NULL}.
#'
#' @param prior_inclusions_prob An (p by 1) vector of prior inclusion
#' probabilities.
#'
#' @param prior_varB A scalar given the (initial) prior variance of the coefficient.
#'
#' @param optimize_prior_varB The optimization method to use for fitting the prior
#' variance.
#'
#' @param check_null_threshold Scalar specifying the threshold on the log-scale
#' to compare likelihood between current estimate and zero (null).
#'
#' @param mle_var_estimator The estimation method of the variance of the MLEs.
#' \code{"naive"} calculates the variances from the weighted least square model,
#' treating weights as constants.
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
#' \item{prior_varB}{Prior variance (after optimization if \code{optimize_prior_varB != "none"}).}
#'
#' @importFrom stats dnorm
#' @importFrom stats uniroot
#' @importFrom stats optim
#' @importFrom Matrix colSums
#' @importFrom sandwich vcovHC
#'
#' @keywords internal
#'
weighted_single_effect_regression <-
  function(y, X,
           W = NULL,
           prior_varB,
           residual_variance = 1,
           prior_inclusion_prob = NULL,
           optimize_prior_varB = c("none", "optim", "EM", "simple"),
           check_null_threshold = 0) {

    optimize_prior_varB <- match.arg(optimize_prior_varB)

    p <- ncol(X)
    n <- nrow(X)

    # Check W
    if (is.null(W)) {
      W <- rep(1, times = n)
    } else if (length(W) == n) {  ## (n by 1) vector
      W <- as.numeric(W)
    } else {
      stop("Dimensions of W and X do not match!")
    }

    ## Assume each column of X has been centralized...
    # Scale X: centralize X by its weighted mean.
    # should I compute the weighted means here? Or in the update_each_effect()?
    X_ <- sweep(X, 2, colWeightedMeans(x = X, w = W), "-")
    ## How about standardization?
    # X_ <- sweep(X_, 2, attr(X, "scaled:scale"), "/")

    # Scale y: y is a n-dim matrix.
    y_ <- y - weighted.mean(y, w = W)

    wxy <- colSums(sweep(X_, 1, W * y_, "*"))
    wx2 <- colSums(sweep(X_^2, 1, W, "*"))
    betahat <- wxy / wx2
    shat2 <- residual_variance / wx2
    # shat2 <- 1/wx2

    # Check prior_inclusion_probability
    ## If `prior_inclusion_probability` is consistent to all variables,
    ## its value should not affect the posterior inclusion probability (alpha).
    if (is.null(prior_inclusion_prob)) {
      prior_inclusion_prob <- rep(1/p, times = p)
    } else {
      if (length(prior_inclusion_prob) != p)
        stop("The length of prior inclusion probability does not match")
    }

    if (optimize_prior_varB != "EM" && optimize_prior_varB != "none") {
      prior_varB <- optimize_prior_variance(optimize_prior_varB, betahat, shat2,
                                  prior_inclusion_prob,
                                  alpha = NULL, mu2 = NULL,
                                  prior_varB_init = prior_varB,
                                  check_null_threshold = check_null_threshold)

      prior_varB <- pmax(prior_varB, 1e-10)
    }

    # compute log of Bayes factor (logBF) and log of posterior odds (logPO)
    # .loge <- function(t) log(t+.Machine$double.eps) (in `elbo.R`)
    zscore2 <- betahat^2 / shat2
    logBF <- 1/2 * (.loge(shat2 /(shat2 + prior_varB)) +
                      zscore2 * prior_varB / (prior_varB + shat2))
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
    post_varB <- (1/shat2 + 1/prior_varB)^(-1)        # posterior variance
    mu <- prior_varB / (prior_varB + shat2) * betahat  # posterior mean
    mu2 <- post_varB + mu^2                      # posterior second moment

    # ABF for WSER model
    logBF_model <- maxlogPO + log(sum(logPO_weighted))
    # = log(sum(ABF x prior_weights))

    if (optimize_prior_varB == "EM") {
      prior_varB <- optimize_prior_variance(optimize_prior_varB, betahat, shat2,
                                   prior_inclusion_prob,
                                   alpha, mu2,
                                   check_null_threshold = check_null_threshold)
      prior_varB <- pmax(prior_varB, 1e-10)
    }
    print(prior_varB)

    return(list(alpha = alpha,      # posterior inclusion probability
                mu = mu,            # posterior mean
                mu2 = mu2,          # posterior second moment
                betahat = betahat,  # maximum likelihood estimator (MLE)
                logBF = logBF,      # layer-wise log of Bayes factor
                logBF_model = logBF_model,  # log of Bayes factor of model
                prior_varB = prior_varB,  # prior variance of coefficients B
    ))
  }


