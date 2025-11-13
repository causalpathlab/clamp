#' Estimate prior variance
#'
#' In this function, deltahat represents the MLE,
#' and shat2 represents the corresponding variance.
#'
#' @importFrom stats optim
#'
#' @keywords internal
optimize_prior_variance <- function(optimize_prior_varD, deltahat, shat2,
                                    prior_inclusion_prob,
                                    alpha = NULL, mu2 = NULL,
                                    prior_varD_init = NULL,
                                    check_null_threshold = 0) {
  prior_varD = prior_varD_init

  if (optimize_prior_varD != "simple") {
    if (optimize_prior_varD == "optim") {
      log_prior_varD <- optim(par = log(max(c(deltahat^2 - shat2, 1), na.rm = T)),
                              fn = neg.optimfunc.logscale,
                              deltahat = deltahat, shat2 = shat2,
                              prior_inclusion_prob = prior_inclusion_prob,
                              method = "Brent", lower = -30, upper = 15)$par
      # If the estimated one is worse than the current one, don't change it
      if (neg.optimfunc.logscale(log_prior_varD, deltahat = deltahat, shat2 = shat2,
                                 prior_inclusion_prob = prior_inclusion_prob) >
          neg.optimfunc.logscale(log(prior_varD), deltahat = deltahat, shat2 = shat2,
                                 prior_inclusion_prob = prior_inclusion_prob)) {
        log_prior_varD <- log(prior_varD)
      }
      prior_varD <- exp(log_prior_varD)
    }
    else if (optimize_prior_varD == "EM") {
      prior_varD <- sum(alpha * mu2)  # sum of second-order of beta_js
    } else
      stop("Invalid option for optimize_prior_varD method")
  }
  ## if (optimize_prior_varD == "simple"), prior_varD is compared with `check_null_threshold`
  ## without any other updates;
  ## if (optimize_prior_varD == "none"), prior_varD is always the pre-assigned coefficient prior
  ## variance and is not compared with any other values.

  ## Set prior_varD exactly 0 if that beats the numerical value by check_null_threshold
  ## in loglik. It means that for parsimony reasons we set estimate of log_prior_varD to zero
  ## if its numerical estimate is only "negligibly" different from zero.
  if (optimfunc.logscale(0, deltahat, shat2, prior_inclusion_prob) +
      check_null_threshold >=
      optimfunc.logscale(prior_varD, deltahat, shat2, prior_inclusion_prob))
    prior_varD <- 0

  return(prior_varD)
}


#' The following is the log-scale of the optimization goal
#' as a function of prior variance V.
#' @keywords internal
optimfunc.logscale <- function(prior_varD, deltahat, shat2,
                               prior_inclusion_prob) {

  # compute log-ABF for each variable
  zscore2 <- deltahat^2 / shat2
  logBF <- 1/2 * (.loge(shat2 /(shat2 + prior_varD)) +
                    zscore2 * prior_varD / (prior_varD + shat2))
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
neg.optimfunc.logscale <- function(log_prior_varD, deltahat,
                                   shat2, prior_inclusion_prob) {
  return(-optimfunc.logscale(exp(log_prior_varD), deltahat, shat2,
                             prior_inclusion_prob))
}
