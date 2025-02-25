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
                                    alpha = NULL, mu2 = NULL,
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
      prior_varB <- sum(alpha * mu2)  # sum of second-order of beta_js
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
  logBF <- 1/2 * (.loge(shat2 /(shat2 + prior_varB)) +
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
