#' @title Update each effect once in linear model; each sub-model is a weighted SER.
#' @param X An (n by p) matrix of regressor variables
#' @param y An n vector of response variable
#' @param s A clamp fit
#' @param W An (n by p) matrix of weights. If \code{W=NULL}, it reduces to an SER
#' @param estimate_prior_variance boolean indicating whether to
#'   estimate prior variance
#' @param check_null_threshold Float, a threshold on the log scale to
#'   compare likelihood between current estimate and zero the null
#' @param robust_method
#' @param robust_estimator
#'
update_each_effect <- function (X, y, s, W=NULL,
                                estimate_prior_variance = FALSE,
                                estimate_prior_method = "optim",
                                check_null_threshold,
                                robust_method = c("none", "huber"),
                                robust_estimator = c("M", "S")) {

  if (!estimate_prior_variance) estimate_prior_method = "none"

  robust_method    <- match.arg(robust_method)
  robust_estimator <- match.arg(robust_estimator)

  # Check weight matrix W
  if (!is.null(W) &
      !(length(W) %in% c(nrow(X), nrow(X) * ncol(X))))
    stop("The dimension of W does not match with that of input X!")

  # Robust estimation regarding W/residuals.
  # The importance weights are assigned to each observation.
  # `s$importance_weight` should be an (n by 1) vector.

  current_R <- y - s$Xr

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
  if (!is.null(W)) {
    W <- sweep(W, 1, s$importance_weight, "*")
  } else {  ## is.null(W)
    W <- s$importance_weight
  }


  # Repeat for each effect to update.
  maxL = nrow(s$alpha)
  if (maxL > 0)
    for (l in 1:maxL) {

      # Residuals belonging to layer l
      Rl <- current_R + compute_Xb(X, s$alpha[l,] * s$mu[l,])

      # Fit WSER
      res <- weighted_single_effect_regression(y=Rl, X=X,
                                  W = W,
                                  residual_variance = s$sigma2,
                                  prior_inclusion_prob = s$pie,
                                  prior_varB = s$prior_varB[l],
                                  optimize_prior_varB = estimate_prior_method,
                                  check_null_threshold = check_null_threshold)

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
      current_R <- Rl - compute_Xb(X,s$alpha[l,] * s$mu[l,])
    }

  s$Xr <- compute_Xb(X, colSums(s$alpha * s$mu))  # linear predictor

  return(s)
}
