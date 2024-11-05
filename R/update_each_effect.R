#' @title Update each effect once in a linear model; each sub-model is a weighted SER.
#'
#' @param X An (n by p) matrix of regressor variables
#'
#' @param y An n vector of response variable
#'
#' @param s A clamp fit
#'
#' @param W An (n by p) matrix of weights. If \code{W=NULL}, it reduces to an SER
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
#' @param estimate_prior_variance boolean indicating whether to
#'   estimate prior variance
#'
#' @param check_null_threshold Float, a threshold on the log scale to
#'   compare likelihood between current estimate and zero the null
#'
#' @param abnormal_proportion a value between 0 and 1. If the number of detected
#'   abnormal subjects exceeds \eqn{abnormal_proportion * nrow(X)},
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
update_each_effect <- function (X, y, s, W=NULL, model,
                  mle_variance_estimator = c("bootstrap", "sandwich", "naive"),
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
  mle_variance_estimator <- match.arg(mle_variance_estimator)

  # Check weight matrix W
  if (!is.null(W) &
      !( (length(W) == nrow(X)) | all(dim(W) == dim(X))  ) )
    stop("The dimension of W does not match with that of input X!")

  # Compute the residuals
  if (s$family == "linear") {

    ## compute the residuals
    current_R <- y - s$Xr

    ## IRLS is not applicable for linear models,
    ## for the consistency of codes, let irls_weight = 1
    irls_weight <- rep(1, times = length(y))

  } else { ## if (s$family %in% c("logistic", "poisson"))

    # Iterative reweighted least-squared (IRLS) for generalized linear models

    ## update the pseudo-response
    psd_rsp <- model$pseudo_response(s$Xr, y)

    ## update the overall log-pseudo-variance
    log_psd_var <- model$log_pseudo_var(s$Xr)
    ## a n-dim vector of weights yielded from IRLS
    irls_weight <- exp(-log_psd_var)
    ## check if there are any abnormal points based on irls_weight
    s$abnormal_subjects <- check_abnormal_subjects(irls_weight)

    ## compute the residuals
    current_R <- psd_rsp - s$Xr
  }

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
    WW <- sweep(W, 1, irls_weight * s$importance_weight, "*")
  } else if (is.vector(W)) {
    WW <- W * irls_weight * s$importance_weight
  } else {  ## is.null(W)
    WW <- irls_weight * s$importance_weight
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
      WW_sub         <- remove_abnormal_subjects(s$abnormal_subjects, WW)

      ## inherit the attributes of X_sub from X
      attr(X_sub, "scaled:scale")  <- attr(X, "scaled:scale")
      attr(X_sub, "scaled:center") <- attr(X, "scaled:center")
    }


    # Update each effect
    for (l in 1:maxL) {

      ## residuals belonging to layer l
      Rl <- current_R + compute_Xb(X_sub, s$alpha[l,] * s$mu[l,])

      ## fit the WSER
      res <- weighted_single_effect_regression(y=Rl, X=X_sub,
                                W = WW_sub,
                                residual_variance = s$sigma2,
                                prior_inclusion_prob = s$pie,
                                prior_varB = s$prior_varB[l],
                                optimize_prior_varB = estimate_prior_method,
                                check_null_threshold = check_null_threshold,
                                mle_variance_estimator = mle_variance_estimator,
                                nboots = nboots,
                                seed = seed)


      # Update the variational estimate of the posterior distributions.
      s$mu[l,]        = res$mu
      s$mu2[l,]       = res$mu2
      s$alpha[l,]     = res$alpha
      s$betahat[l,]   = res$betahat
      s$prior_varB[l] = res$prior_varB
      s$logBF[l]      = res$logBF_model
      s$logBF_variable[l,] = res$logBF
      # s$KL[l]     = -res$loglik +
      #   SER_posterior_e_loglik(X,R,s$sigma2,res$alpha * res$mu,
      #                          res$alpha * res$mu2)

      # Update the current residuals
      current_R <- Rl - compute_Xb(X_sub, s$alpha[l,] * s$mu[l,])
    }

    s$Xr <- compute_Xb(X, colSums(s$alpha * s$mu))  # update linear predictor
  }

  return(s)
}
