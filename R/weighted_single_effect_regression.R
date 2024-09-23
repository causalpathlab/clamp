#' @rdname weighted_single_effect_regresion
#'
#' @title Weighted single effect regression
#'
#' @description
#' The WSER function is to compute the posterior distribution of the regression
#' coefficients of a WSER model. Reference:
#' <https://github.com/stephenslab/susieR/blob/master/R/single_effect_regression.R>
#'
#' @param y an (n by 1) vector of responses.
#'
#' @param X An (n by p) data matrix.
#'
#' @param W An (n by 1) vector or an (n by p) matrix.
#' If \code{W} is an (n by 1) vector, it contains the W of each subject;
#' for a generalized linear model using iterative reweighted least squared
#' approach, \eqn{W = exp(-logw2)};
#' If \code{W} is an (n by p) matrix, each entry is the weight of each each
#' entry of \eqn{X}; it should be of the same size as \eqn{X}.
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

    # Check W
    if (is.null(W)) {
      W <- rep(1, times = nrow(X))
    } else if (length(W) == nrow(X)) {  ## (n by 1) vector
      W <- as.numeric(W)
    } else if (length(W) == nrow(X)*ncol(X)) {  ## (n by p) matrix
      if (dim(W)[1] != nrow(X))
        stop("Dimensions of the W do not match!")
    } else {
      stop("Dimensions of the W do not match!")
    }

    # Scale X
    X_ <- sweep(X, 2, attr(X, "scaled:center"), "-")  ## centralized
    X_ <- sweep(X_, 2, attr(X, "scaled:scale"), "/")

    # Update the MLE (using scaled X)
    if (length(W) == nrow(X)) {  ## weight is an (n by 1) vector
      XtWX <- colSums(sweep(X_ * X_, 1, W, "*"))
      XtWy <- colSums(sweep(X_, 1, W * y, "*"))
    } else {  ## weight is an (n by p) matrix
      XtWX <- colSums(X_ * X_ * W)
      XtWy <- colSums(sweep(X_ * W, 1, y, "*"))
    }
    ## [Wrap up these matrix computation.]!

    shat2 <- residual_variance / XtWX
    betahat <- XtWy / XtWX  # = shat2 * XtWy

    # Check prior_inclusion_probability
    ## If `prior_inclusion_probability` is consistent to all variables,
    ## its value should not affect the posterior inclusion probability (alpha).
    if (is.null(prior_inclusion_prob)) {
      prior_inclusion_prob <- rep(1/p, times = p)  ## other value?
    } else {
      if (length(prior_inclusion_prob) != p)
        stop("The length of prior inclusion probability does not match")
    }


    if (optimize_prior_varB != "EM" && optimize_prior_varB != "none") {
      prior_varB <- optimize_prior_variance(optimize_prior_varB, betahat, shat2,
                                  prior_inclusion_prob,
                                  alpha = NULL, post_mean2 = NULL,
                                  prior_varB_init = prior_varB,
                                  check_null_threshold = check_null_threshold)
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
    post_var <- 1 / (1/shat2 + 1/prior_varB)                      # posterior variance
    post_mean <- prior_varB / (prior_varB + shat2) * XtWy / XtWX  # posterior mean
    post_mean2 <- post_var + post_mean^2                          # posterior second moment

    # ABF for WSER model
    logBF_model <- maxlogPO + log(sum(logPO_weighted))
    # = log(sum(ABF x prior_weights))

    if (optimize_prior_varB == "EM") {
      prior_varB <- optimize_prior_variance(optimize_prior_varB, betahat, shat2,
                                   prior_inclusion_prob,
                                   alpha, post_mean2,
                                   check_null_threshold = check_null_threshold)
      }

    return(list(alpha = alpha,      # posterior inclusion probability
                mu = post_mean,     # posterior mean
                mu2 = post_mean2,   # posterior second moment
                betahat = betahat,  # maximum likelihood estimator (MLE)
                # shat2 = shat2,      # variance of MLE
                logBF = logBF,      # layer-wise log of Bayes factor
                logBF_model = logBF_model,  # log of Bayes factor of model
                prior_varB = prior_varB  # prior variance of coefficients B
    ))
  }


#' Estimate prior variance
#'
#' In this function, betahat represents the MLE,
#' and shat2 represents the corresponding variance.
#'
#' @importFrom stats optim
#'
#' @keywords internal
optimize_prior_variance <- function(optimize_prior_varB, betahat, shat2,
                                    prior_inclusion_prob,
                                    alpha = NULL, post_mean2 = NULL,
                                    prior_varB_init = NULL,
                                    check_null_threshold = 0) {
  prior_varB = prior_varB_init

  if (optimize_prior_varB != "simple") {
    if (optimize_prior_varB == "optim") {
      log_prior_varB <- optim(par = log(max(c(betahat^2 - shat2, 1), na.rm = T)),
                    fn = neg.optimfunc.logscale,
                    betahat = betahat, shat2 = shat2,
                    prior_inclusion_prob = prior_inclusion_prob,
                    method = "Brent", lower = -30, upper = 15)$par
      # If the estimated one is worse than the current one, don't change it
      if (neg.optimfunc.logscale(log_prior_varB, betahat = betahat, shat2 = shat2,
                                 prior_inclusion_prob = prior_inclusion_prob) >
          neg.optimfunc.logscale(log(prior_varB), betahat = betahat, shat2 = shat2,
                                 prior_inclusion_prob = prior_inclusion_prob)) {
        log_prior_varB <- log(prior_varB)
      }
      prior_varB <- exp(log_prior_varB)
    }
    else if (optimize_prior_varB == "EM") {
      prior_varB <- sum(alpha * post_mean2)  # second-order of beta_js (WHY?!!)
    } else
      stop("Invalid option for optimize_prior_varB method")
  }
  ## if (optimize_prior_varB == "simple"), prior_varB is compared with `check_null_threshold`
  ## without any other updates;
  ## if (optimize_prior_varB == "none"), prior_varB is always the pre-assigned coefficient prior
  ## variance and is not compared with any other values.

  ## Set prior_varB exactly 0 if that beats the numerical value by check_null_threshold
  ## in loglik. It means that for parsimony reasons we set estimate of log_prior_varB to zero
  ## if its numerical estimate is only "negligibly" different from zero.
  if (optimfunc.logscale(0, betahat, shat2, prior_inclusion_prob) +
      check_null_threshold >=
      optimfunc.logscale(prior_varB, betahat, shat2, prior_inclusion_prob))
    prior_varB <- 0

  return(prior_varB)
}


#' The following is the log-scale of the optimization goal
#' as a function of prior variance V.
#' @keywords internal
optimfunc.logscale <- function(prior_varB, betahat, shat2,
                               prior_inclusion_prob) {

  # compute log-ABF for each variable
  zscore2 <- betahat^2 / shat2
  logBF <- 1 / 2 * ((log(shat2) - log(shat2 + prior_varB)) +
                       zscore2 * prior_varB / (prior_varB + shat2))
  # deal with special case of infinite shat2 (e.g. happens if X does not vary)
  logBF[is.infinite(shat2)] <- 0
  maxlogBF <- max(logBF)

  # w is proportional to ABF, but subtract max for numerical stability
  w <- exp(logBF - maxlogBF)

  # Update the posterior estimates
  # Posterior prob for each variable
  w_weighted <- w * prior_inclusion_prob
  sum_w_weighted <- sum(w_weighted)

  obj <- log(sum_w_weighted) + maxlogBF  # logBF for a WSER

  return(obj)
}


#' The following is the negative of the objective function
#' @keywords internal
neg.optimfunc.logscale <- function(log_prior_varB, betahat,
                                   shat2, prior_inclusion_prob) {
  return(-optimfunc.logscale(exp(log_prior_varB), betahat, shat2,
                             prior_inclusion_prob))
}
