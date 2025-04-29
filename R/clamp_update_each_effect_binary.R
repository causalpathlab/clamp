#' @title Update each effect once in a linear model; each sub-model is a
#' IPW-SER.
#'
#' @param X An (n by p) matrix of regressor variables
#'
#' @param y An n vector of response variable
#'
#' @param s A clamp fit
#'
#' @param W An (n by p) matrix of weights. If \code{W=NULL}, it reduces to an
#' SER
#'
#' @param mle_estimator The estimation method of the ATE estimator.
#' \code{"mHT"} applies the modified Horvitz-Thompson estimator, and
#' \code{"WLS"} applies the equivalent weighted least-squares estimator.
#'
#' @param mle_variance_estimator The estimation method of the variance of the
#' MLEs, or equivalently, the variance of (Horvitz-Thompson) ATE estimators.
#' \code{"bootstrap"} estimates the variances from \code{nboots} bootstrap
#' replicates,
#' \code{"sandwich"} uses sandwich robust variance estimators.
#' \code{"naive"} calculates the variances from the weighted least square model,
#' treating weights as constants (not recommended).
#'
#' @param nboots The number of bootstrap replicates. By default,
#' \code{nboots=100}.
#'
#' @param seed Random seed utilized when
#' \code{mle_variance_estimator="bootstrap"}.
#' If \code{seed=NULL}, its default is \code{Sys.time()}.
#'
#' @param standardize logical. If \code{standardize=TRUE}, then each column of
#'. X will be standardized by its corresponding weighted sd.
#'
#' @param estimate_prior_variance boolean indicating whether to
#'   estimate prior variance
#'
#' @param check_null_threshold Float, a threshold on the log scale to
#'   compare likelihood between current estimate and zero the null
#'
#' @param abnormal_proportion a value between 0 and 1. If the number of detected
#'   abnormal subjects exceeds \code{abnormal_proportion * nrow(X)},
#'   stop fitting the model.
#'
#' @param robust_method A string, whether and which robust method is applied
#'   when fitting the model. \code{robust_method="none"} specifies that no
#'   robust method is applied. \code{robust_method="huber"} specifies that the
#'   Huber weighting method is applied.
#'
#' @param robust_estimator A string, which robust estimator is applied.
#'   \code{robust_estimator="M"} indicates the M-estimator is applied, and
#'   \code{robust_estimator="S"} indicates the S-estimator is applied.
#'
clamp_update_each_effect_binary <- function (X, y, s, W=NULL,
                  mle_estimator = c("mHT", "WLS"),
                  mle_variance_estimator = c("bootstrap", "sandwich"),
                  nboots = 100,
                  seed = NULL,
                  estimate_prior_variance = FALSE,
                  estimate_prior_method = "optim",
                  check_null_threshold = 0,
                  abnormal_proportion = 0.5,
                  robust_method = c("none", "huber"),
                  robust_estimator = c("M", "S")) {

  if (!estimate_prior_variance) estimate_prior_method = "none"

  robust_method    <- match.arg(robust_method)
  robust_estimator <- match.arg(robust_estimator)
  mle_estimator    <- match.arg(mle_estimator)
  mle_variance_estimator <- match.arg(mle_variance_estimator)

  # Check weight matrix W
  if (!is.null(W) &
      !( (length(W) == nrow(X)) | all(dim(W) == dim(X))  ) )
    stop("The dimension of W does not match with that of input X!")


  # Compute the residuals
  current_R <- y - s$Xr

  # Robust estimation regarding W/residuals.
  # The importance weights are assigned to each observation.
  # `s$importance_weight` should be an (n by 1) vector.
  if (robust_method != "none" & robust_estimator == "S"){

    # Assign importance weight based on the residuals yielded from the
    # previous iteration.
    s$importance_weight <- robust_importance_weights(current_R,
                            robust_method = robust_method,
                            robust_estimator = robust_estimator,
                            previous_importance_weight = s$importance_weight)
  } else {
    ## robust_method == "none" | robust_estimator=="M"
    s$importance_weight <- robust_importance_weights(current_R,
                              robust_method = robust_method,
                              robust_estimator = robust_estimator)
  }

  # Define the weight (matrix)
  if (is.matrix(W)) {
    WW <- sweep(W, 1, s$importance_weight, "*")
  } else if (is.vector(W)) {
    WW <- W * s$importance_weight
  } else {
    WW <- s$importance_weight
  }


  # Repeat for each effect to update.
  maxL = nrow(s$alpha)
  if (maxL > 0) {

    # Remove abnormal points
    if (length(s$abnormal_subjects) > abnormal_proportion * nrow(X)) {
      stop("Too many abnormal subjected detected...")

    } else {
      X_sub         <- remove_abnormal_subjects(s$abnormal_subjects, X)
      current_R_sub <- remove_abnormal_subjects(s$abnormal_subjects, current_R)
      WW_sub        <- remove_abnormal_subjects(s$abnormal_subjects, WW)

      ## inherit the attributes of X_sub from X
      attr(X_sub, "scaled:scale")  <- attr(X, "scaled:scale")
      attr(X_sub, "scaled:center") <- attr(X, "scaled:center")
    }

    # Update each effect
    for (l in 1:maxL) {

      ## residuals belonging to layer l
      Rl <- current_R + compute_Xb(X = X_sub, b = s$alpha[l,]*s$mu[l,])


      ## fit the ipw-ser
      res <- ipw_single_effect_regression_binary(y=Rl, X=X_sub,
                                W = WW_sub,
                                residual_variance = s$sigma2,
                                prior_inclusion_prob = s$pie,
                                prior_varD = s$prior_varD[l],
                                optimize_prior_varD = estimate_prior_method,
                                check_null_threshold = check_null_threshold,
                                mle_estimator = mle_estimator,
                                mle_variance_estimator = mle_variance_estimator,
                                nboots = nboots,
                                seed = seed)


      # Update the variational estimate of the posterior distributions.
      s$mu[l,]        = res$mu
      s$mu2[l,]       = res$mu2
      s$alpha[l,]     = res$alpha
      s$deltahat[l,]   = res$deltahat
      s$shat2[l,]     = res$shat2
      s$prior_varD[l] = res$prior_varD
      s$logBF[l]      = res$logBF_model
      s$logBF_variable[l,] = res$logBF
      s$EWR2[l]       = res$EWR2
      s$Eloglik[l]    = res$Eloglik

      # Update the current residuals
      current_R <- Rl - compute_Xb(X_sub, s$alpha[l,] * s$mu[l,])
    }

    s$Xr <- compute_Xb(X, colSums(s$alpha * s$mu))  # update linear predictor
  }

  return(s)
}
