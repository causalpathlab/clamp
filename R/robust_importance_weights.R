#' @title Robust estimations
#' Copied from gsusie
#'
#' @description Outliers or some unexpected resid due to extremely small/large
#' weights may occur in both observations and intermediate
#' procedures. Hence, robust estimations are be considered to "down-weigh"
#' the abnormal observations.
#'
#' @param resid a vector that may contains outliers.
#' May be weights (inverse of pseudo-variance) or residuals(?!)
#'
#' @param robust_method If \code{method = "none"}, then all subjects are
#' included in the coefficient estimation.
#' If \code{method = "huber"}, huber weights are
#' assigned to each subject based on the \dQuote{pseudo-residuals} \code{rr}.
#'
#' @param robust_estimator If \code{robust_estimator = "M"},
#' M-estimation is performed. If \code{robust_estimator = "S"},
#' S-estimation is performed. This argument specifies the tuning method,
#' i.e., the method for defining outliers in each iteration, when applying
#' Huber or Bisquare reweighting method (\code{robust_method = "huber"} or
#' \code{robust_method = "bisquare"}).
#'
#' @param previous_importance_weight importance weights in the previous iteration,
#' used to calculate the importance weights in the current iteration for the
#' S-estimation. If \code{previous_importance_weight} is NULL, it indicates the
#' first iteration, and thus S-estimation need to be initialized.
#'
#' @param huber_tuning_k The tuning parameter of the Huber's M-estimator.
#'
#' @returns {importance_weight} An vector of length n, the importance weights.
#' If a subject is to be removed from the current iteration,
#' its importance weights is 0.
#'
#' @importFrom stats quantile
#'
#' @keywords internal
#'

robust_importance_weights <- function(
    resid,
    robust_method = c("none", "huber"),
    robust_estimator = c("M", "S"),
    previous_importance_weight = NULL,
    huber_tuning_k = 1.345)
  {

  robust_method <- match.arg(robust_method)
  robust_estimator <- match.arg(robust_estimator)

  if (robust_method == "none"){

    importance_weight <- rep(1, times = length(resid))
  }
  else if (robust_method == "huber") {

    if (robust_estimator == "M"){  # M-Estimation

      hat_sd_resid <- median(abs(resid - median(resid))) / 0.6745
      std_resid <- resid /  hat_sd_resid  # standardized residuals u
      importance_weight <- huber_weight_(std_resid, huber_tuning_k)
    }
    else {  # S-Estimation

      if (is.null(previous_importance_weight))
      { # iteration = 1, initialization
        hat_sd_resid <- median(abs(resid - median(resid))) / 0.6745
        std_resid <- resid /  hat_sd_resid  # standardized residuals u
        importance_weight <- huber_weight_(std_resid, huber_tuning_k)
      }
      else { # iteration > 1

        # robust estimation of the standard deviation of the residuals.
        hat_sigma_resid <- sqrt(
          sum(sweep(as.matrix(previous_importance_weight), 1, resid^2, "*")) /
            (length(resid) * 0.199))
        std_resid <- resid / hat_sigma_resid
        importance_weight <-
          huber_loss_(std_resid, huber_tuning_k) / std_resid^2
      }
    }
  }

  return(importance_weight)

}

huber_loss_ <- function(r, k) {
  out <- ifelse(abs(r) <= k, r^2/2, k*abs(r) - k^2/2)
  return(out)
}

huber_weight_ <- function(r, k) {
  out <- ifelse(abs(r) <= k, 1, k / abs(r))
  return(out)
}

