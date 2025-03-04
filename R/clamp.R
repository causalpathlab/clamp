#' ...
#'
#' @param X An n by p matrix of covariates.
#'
#' @param y The observed responses, a vector of length n.
#'
#' @param W The weight matrix of size (n by p) or a n-dim vector.
#'   By default, \code{W=NULL}.
#'
#' @param maxL Maximum number of non-zero effects in the susie
#'   regression model. If L is larger than the number of covariates, p,
#'   L is set to p.
#'
#' @param family A description of error distribution and link function used in
#'   the model. \code{family="linear"} stands for a linear regression model
#'   with identity link, \code{family="binomial"} stands for a logistic
#'   regression model, and \code{family="poisson"} stands for a Poisson
#'   regression model.
#'
#' @param scaled_prior_variance The prior variance, divided by
#'   \code{var(y)} (or by \code{(1/(n-1))yty} for
#'   \code{susie_suff_stat}); that is, the prior variance of each
#'   non-zero element of b is \code{var(y) * scaled_prior_variance}. The
#'   value provided should be either a scalar or a vector of length
#'   \code{L}. If \code{estimate_prior_variance = TRUE}, this provides
#'   initial estimates of the prior variances.
#'
#' @param residual_variance Variance of the residual. If
#'   \code{estimate_residual_variance = TRUE}, this value provides the
#'   initial estimate of the residual variance. By default, it is set to
#'   \code{var(y)} in \code{clamp} and \code{(1/(n-1))yty} in
#'   \code{susie_suff_stat}.
#'
#' @param prior_inclusion_prob A vector of length p, in which each entry
#'   gives the prior probability that corresponding column of X has a
#'   nonzero effect on the outcome, y.
#'
#' @param null_weight Prior probability of no effect (a number between
#'   0 and 1, and cannot be exactly 1).
#'
#' @param standardize If \code{standardize = TRUE}, standardize the
#'   columns of X to unit variance prior to fitting (or equivalently
#'   standardize XtX and Xty to have the same effect). Note that
#'   \code{scaled_prior_variance} specifies the prior on the
#'   coefficients of X \emph{after} standardization (if it is
#'   performed). If you do not standardize, you may need to think more
#'   carefully about specifying \code{scaled_prior_variance}. Whatever
#'   your choice, the coefficients returned by \code{coef} are given for
#'   \code{X} on the original input scale. Any column of \code{X} that
#'   has zero variance is not standardized.
#'
#' @param intercept If \code{intercept = TRUE}, the intercept is
#'   fitted; it \code{intercept = FALSE}, the intercept is set to
#'   zero. Setting \code{intercept = FALSE} is generally not
#'   recommended.
#'
#' @param estimate_residual_variance If
#'   \code{estimate_residual_variance = TRUE}, the residual variance is
#'   estimated, using \code{residual_variance} as an initial value. If
#'   \code{estimate_residual_variance = FALSE}, the residual variance is
#'   fixed to the value supplied by \code{residual_variance}.
#'
#' @param estimate_prior_variance If \code{estimate_prior_variance =
#'   TRUE}, the prior variance is estimated (this is a separate
#'   parameter for each of the L effects). If provided,
#'   \code{scaled_prior_variance} is then used as an initial value for
#'   the optimization. When \code{estimate_prior_variance = FALSE}, the
#'   prior variance for each of the L effects is determined by the
#'   value supplied to \code{scaled_prior_variance}.
#'
#' @param estimate_prior_method The method used for estimating prior
#'   variance. When \code{estimate_prior_method = "simple"} is used, the
#'   likelihood at the specified prior variance is compared to the
#'   likelihood at a variance of zero, and the setting with the larger
#'   likelihood is retained.
#'
#' @param check_null_threshold When the prior variance is estimated,
#'   compare the estimate with the null, and set the prior variance to
#'   zero unless the log-likelihood using the estimate is larger by this
#'   threshold amount. For example, if you set
#'   \code{check_null_threshold = 0.1}, this will "nudge" the estimate
#'   towards zero when the difference in log-likelihoods is small. A
#'   note of caution that setting this to a value greater than zero may
#'   lead the IBSS fitting procedure to occasionally decrease the ELBO.
#'
#' @param prior_tol When the prior variance is estimated, compare the
#'   estimated value to \code{prior_tol} at the end of the computation,
#'   and exclude a single effect from PIP computation if the estimated
#'   prior variance is smaller than this tolerance value.
#'
#' @param residual_variance_upperbound Upper limit on the estimated
#'   residual variance. It is only relevant when
#'   \code{estimate_residual_variance = TRUE}.
#'
#' @param robust_method
#'
#' @param robust_estimator
#'
#' @param coverage A number between 0 and 1 specifying the
#'   \dQuote{coverage} of the estimated confidence sets.
#'
#' @param min_abs_corr Minimum absolute correlation allowed in a
#'   credible set. The default, 0.5, corresponds to a squared
#'   correlation of 0.25, which is a commonly used threshold for
#'   genotype data in genetic studies.
#'
#' @param na.rm Drop any missing values in y from both X and y.
#'
#' @param max_iter Maximum number of IBSS iterations to perform.
#'
#' @param tol A small, non-negative number specifying the convergence
#'   tolerance for the IBSS fitting procedure. The fitting procedure
#'   will halt when the difference in the variational lower bound, or
#'   \dQuote{ELBO} (the objective function to be maximized), is
#'   less than \code{tol}.
#'
#' @param verbose If \code{verbose = TRUE}, the algorithm's progress,
#'   and a summary of the optimization settings, are printed to the
#'   console.
#'
#' @param track_fit If \code{track_fit = TRUE}, \code{trace}
#'   is also returned containing detailed information about the
#'   estimates at each iteration of the IBSS fitting procedure.
#'
#' @param residual_variance_lowerbound Lower limit on the estimated
#'   residual variance. It is only relevant when
#'   \code{estimate_residual_variance = TRUE}.
#'
#' @param n_purity Passed as argument \code{n_purity} to
#'   \code{\link{clamp_get_cs}}.
#'
#' @return A \code{"clamp"} object with some or all of the following
#'   elements:
#'
#' \item{alpha}{An maxL by p matrix of posterior inclusion probabilites.}
#'
#' \item{mu}{An maxL by p matrix of posterior means, conditional on
#'   inclusion.}
#'
#' \item{mu2}{An maxL by p matrix of posterior second moments,
#'   conditional on inclusion.}
#'
#' \item{betahat}{An maxL by p matrix of maximum likelihood estimator,
#'   conditional on inclusion.}
#'
#' \item{Xr}{A vector of length n, equal to \code{X \%*\% colSums(alpha
#'   * mu)}.}
#'
#' \item{logBF}{log-Bayes Factor for each single effect.}
#'
#' \item{logBF_variable}{log-Bayes Factor for each variable and single effect.}
#'
#' \item{intercept}{Intercept (fixed or estimated).}
#'
#' \item{sigma2}{Residual variance (fixed or estimated).}
#'
#' \item{prior_varB}{Prior variance of the non-zero elements of b, equal to
#'   \code{scaled_prior_variance * var(y)}.}
#'
#' \item{elbo}{The value of the variational lower bound, or
#'   \dQuote{ELBO} (objective function to be maximized), achieved at
#'   each iteration of the IBSS fitting procedure.}
#'
#' \item{importance_weight}{A vector of length n containing the importance
#'   weights computed when applying the robust methods (or all 1 if
#'   \code{robust_method="none"}.)}
#'
#' \item{fitted}{A vector of length n containing the fitted values of
#'   the outcome.}
#'
#' \item{sets}{Credible sets estimated from model fit; see
#'   \code{\link{clamp_get_cs}} for details.}
#'
#' \item{pip}{A vector of length p giving the (marginal) posterior
#'   inclusion probabilities for all p covariates.}
#'
#' \item{z}{A vector of univariate z-scores.}
#'
#' \item{niter}{Number of IBSS iterations that were performed.}
#'
#' \item{converged}{\code{TRUE} or \code{FALSE} indicating whether
#'   the IBSS converged to a solution within the chosen tolerance
#'   level.}

#'
#' @importFrom stats var
#' @importFrom utils modifyList
#'
#' @export
#'
clamp <- function (X, y,
                   W = NULL, ## IPW matrix, should be of same size of X
                   maxL = min(10,ncol(X)),
                   family = c("linear", "binomial", "poisson"),
                   scaled_prior_variance = 0.2,
                   residual_variance = NULL,
                   prior_inclusion_prob = NULL,
                   null_weight = 0,
                   standardize = TRUE,
                   intercept = TRUE,
                   estimate_residual_variance = TRUE,
                   estimate_prior_variance = TRUE,
                   estimate_prior_method = c("optim", "EM", "simple"),
                   check_null_threshold = 0,
                   prior_tol = 1e-9,
                   residual_variance_upperbound = Inf,
                   robust_method = c("none", "huber"),
                   robust_estimator = c("M", "S"),
                   coverage = 0.95,
                   min_abs_corr = 0.5,
                   na.rm = FALSE,
                   max_iter = 500,
                   tol = 1e-3,
                   verbose = FALSE,
                   track_fit = FALSE,
                   residual_variance_lowerbound = var(drop(y))/1e4,
                   n_purity = 100) {

  # Process input estimate_prior_method.
  estimate_prior_method = match.arg(estimate_prior_method)
  # Process input robust_method
  robust_method = match.arg(robust_method)
  # Process input robust_estimator
  robust_estimator = match.arg(robust_estimator)

  # Check input X.
  if (!(is.double(X) & is.matrix(X)) & !inherits(X,"CsparseMatrix") &
      is.null(attr(X,"matrix.type")))
    stop("Input X must be a double-precision matrix, or a sparse matrix, or ",
         "a trend filtering matrix")
  if (is.numeric(null_weight) && null_weight == 0)
    null_weight = NULL
  if (!is.null(null_weight) && is.null(attr(X,"matrix.type"))) {
    if (!is.numeric(null_weight))
      stop("Null weight must be numeric")
    if (null_weight < 0 || null_weight >= 1)
      stop("Null weight must be between 0 and 1")
    if (missing(prior_inclusion_prob))
      prior_inclusion_prob = c(rep(1/ncol(X) * (1 - null_weight),ncol(X)),
                               null_weight)
    else
      prior_inclusion_prob = c(prior_inclusion_prob * (1-null_weight),
                               null_weight)
    X = cbind(X,0)
  }
  if (anyNA(X))
    stop("Input X must not contain missing values")
  if (anyNA(y)) {
    if (na.rm) {
      samples_kept = which(!is.na(y))
      y = y[samples_kept]
      X = X[samples_kept,]
    } else
      stop("Input y must not contain missing values")
  }
  p = ncol(X)
  if (p > 1000 & !requireNamespace("Rfast",quietly = TRUE))
    warning_message("For an X with many columns, please consider installing",
                    "the Rfast package for more efficient credible set (CS)",
                    "calculations.", style='hint')

  # Process input family
  family = match.arg(family)
  model <- list()
  model$type <- family
  model$loglik <- Eloglik_linear
  # switch(family,
  #        "linear" =
  #          {
  #            model$loglik <- Eloglik_linear
  #          },
  #        "binomial" =
  #          {
  #            model$log_psd_var <-  ## logw2
  #            model$psd_resp ## zz
  #            model$loglik <-
  #          },
  #         "poisson" ={}
  #        )

  # Check input y.
  n = nrow(X)
  mean_y = mean(y)

  # Check weight matrix W
  if (!is.null(W) &
      !(length(W) %in% c(nrow(X), nrow(X) * ncol(X))))
    stop("The dim of Weights does not match with the dim of input X!")

  # Center and scale input.
  if (intercept)
    y_ = y - mean_y

  # Set three attributes for matrix X: attr(X,'scaled:center') is a
  # p-vector of column means of X if center=TRUE, a p vector of zeros
  # otherwise; 'attr(X,'scaled:scale') is a p-vector of column
  # standard deviations of X if scale=TRUE, a p vector of ones
  # otherwise; 'attr(X,'d') is a p-vector of column sums of
  # X.standardized^2,' where X.standardized is the matrix X centered
  # by attr(X,'scaled:center') and scaled by attr(X,'scaled:scale').
  out = compute_colstats(X,center = intercept, scale = standardize)
  attr(X,"scaled:center") = out$scaled_center
  attr(X,"scaled:scale")  = out$scaled_scale
  attr(X,"d") = out$d  ## applied in computing ERSS (`get_ER2()`)

  # Initialize susie fit.
  s <- init_setup(n, p, maxL, family,
                  scaled_prior_variance,
                  residual_variance,
                  prior_inclusion_prob,
                  null_weight,
                  as.numeric(var(y)),
                  standardize)
  s <- init_finalize(s)

  # Initialize elbo to NA.
  elbo = rep(as.numeric(NA),max_iter + 1)
  elbo[1] = -Inf;
  tracking = list()
  # Initialize log-likelihood into NA
  loglik = rep(as.numeric(NA),max_iter+1)  ## not required to know...
  loglik[1] = -Inf

  for (tt in 1:max_iter) {

    if (track_fit)
      tracking[[tt]] = clamp_slim(s)

    s <- update_each_effect(X=X, y=y_, s=s, W=W,
                            estimate_prior_variance = estimate_prior_variance,
                            estimate_prior_method = estimate_prior_method,
                            check_null_threshold = check_null_threshold,
                            robust_method = robust_method,
                            robust_estimator = robust_estimator)

    # Compute objective before updating residual variance because part
    # of the objective s$kl has already been computed under the
    # residual variance before the update.
    elbo[tt+1] = get_elbo(X, y_, s, model=model)
    loglik[tt+1] = model$loglik(X, y_, s)

    if (verbose){
      print(paste0("#iteration:", tt))
      print(paste0("objective:", elbo[tt+1]))
    }

    # if (abs(elbo[tt+1] - elbo[tt]) < tol) {
    if ((elbo[tt+1] - elbo[tt]) < tol) {
      s$converged = TRUE
      break
    }

    if (estimate_residual_variance) {
      s$sigma2 = pmax(residual_variance_lowerbound,
                      estimate_residual_variance_fun(X,y_,s))
      if (s$sigma2 > residual_variance_upperbound)
        s$sigma2 = residual_variance_upperbound
    }

  }

  # Remove first (infinite) entry, and trailing NAs.
  elbo = elbo[2:(tt+1)]
  s$elbo = elbo
  loglik = loglik[2:(tt+1)]
  s$loglik = loglik

  s$niter = tt

  if (is.null(s$converged)) {
    warning(paste("IBSS algorithm did not converge in",max_iter,"iterations!"))
    s$converged = FALSE
  }

  if (intercept) {  ################### How about the glm cases?
    # Estimate unshrunk intercept.
    s$intercept = mean_y - sum(attr(X,"scaled:center") *
                        (colSums(s$alpha * s$mu)/attr(X,"scaled:scale")))
    s$fitted = s$Xr + mean_y
  } else {
    s$intercept = 0
    s$fitted = s$Xr
  }
  s$fitted = drop(s$fitted)
  names(s$fitted) = `if`(is.null(names(y)),rownames(X),names(y))

  if (track_fit)
    s$trace = tracking

  # SuSiE CS and PIP.
  if (!is.null(coverage) && !is.null(min_abs_corr)) {
    s$sets = clamp_get_cs(s,coverage = coverage,X = X,
                                  min_abs_corr = min_abs_corr,
                                  n_purity = n_purity)
    s$pip = clamp_get_pip(s,prune_by_cs = FALSE,prior_tol = prior_tol)
  }

  if (!is.null(colnames(X))) {
    variable_names = colnames(X)
    if (!is.null(null_weight)) {
      variable_names[length(variable_names)] = "null"
      names(s$pip) = variable_names[-p]
    } else
      names(s$pip)    = variable_names
    colnames(s$alpha) = variable_names
    colnames(s$mu)    = variable_names
    colnames(s$mu2)   = variable_names
    colnames(s$betahat) = variable_names
    colnames(s$logBF_variable) = variable_names
  }

  # For prediction.
  s$X_column_scale_factors = attr(X,"scaled:scale")

  return(s)
}
