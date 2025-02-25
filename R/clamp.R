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
#'   with identity link, \code{family="logistic"} stands for a logistic
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
#' @param mle_variance_estimator The method to estimate the variances of
#'   IPW estimators. \code{"bootstrap"} estimates the variance via 100
#'   bootstrap replicates, \code{"sandwich"} computes the robust sandwich
#'   variance estimates, and \code{"naive"} calculates the variances based
#'   on the corresponding weighted linear regression models.
#'
#' @param burn_in A number, referring to the number of initial iterations
#'   being discarded. It is used to check the convergence when
#'   \code{mle_variance_estimator="bootstrap"}. Some pilot simulations are
#'   recommended to determine the burn-in period.
#'
#' @param converge_window A number. If, after the burn-in period, the
#'   average of ELBO of \code{converge_window} successive iterations converges,
#'   then the algorithm is considered converged.
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
#' @param abnormal_proportion a value between 0 and 1. If the number of detected
#'   abnormal subjects exceeds \eqn{abnormal_proportion * nrow(X)},
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
#' @importFrom stats weighted.mean
#' @importFrom utils modifyList
#'
#' @export
#'
clamp <- function (X, y,
                   W = NULL, ## IPW matrix, should be of same size of X
                   maxL = min(10,ncol(X)),
                   # family = "linear",
                   scaled_prior_variance = 0.2,
                   residual_variance = NULL,
                   prior_inclusion_prob = NULL,
                   standardize = FALSE,
                   intercept = TRUE,
                   estimate_residual_variance = TRUE,
                   estimate_prior_variance = TRUE,
                   estimate_prior_method = c("optim", "EM", "simple"),
                   check_null_threshold = 0,
                   mle_estimator = c("mHT", "WLS"),
                   mle_variance_estimator = c("bootstrap", "sandwich", "naive"),
                   nboots = 100,
                   seed = NULL,
                   burn_in = 10,
                   converge_window = 10,
                   prior_tol = 1e-9,
                   residual_variance_upperbound = Inf,
                   abnormal_proportion = 0.5,
                   robust_method = c("none", "huber"),
                   robust_estimator = c("M", "S"),
                   coverage = 0.95,
                   min_abs_corr = 0.5,
                   na.rm = FALSE,
                   max_iter = 200,
                   tol = 1e-2,
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

  # Process input mle_variance_estimator
  mle_variance_estimator = match.arg(mle_variance_estimator)

  # # Process input family and the related functions
  # family = match.arg(family)
  # model <- list()
  # model$family <- family

  # Check input X.
  if (!(is.double(X) & is.matrix(X)) & !inherits(X,"CsparseMatrix") &
      is.null(attr(X,"matrix.type")))
    stop("Input X must be a double-precision matrix, or a sparse matrix, or ",
         "a trend filtering matrix")
  ## sparse matrix?!

  # Check input.
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

  ## Check the column names of input X
  if (is.null(colnames(X)))
    colnames(X) <- paste0("X", 1:ncol(X))


  # Center and scale input.
  if (mle_estimator == "mHT") {
    # When applying the modified Horvitz-Thompson estimator,
    # the input X and response Y does not need to be scaled.
    # Each entry of X should be either 0 or 1.
    out = compute_colstats(X, center = F, scale = F)
  } else if (mle_estimator == "WLS") {
    out = compute_colstats(X, W, center = intercept, scale = standardize)
  }
  # Set three attributes for matrix X: attr(X,'scaled:center') is a
  # p-vector of column means of X if center=TRUE, a p vector of zeros
  # otherwise; 'attr(X,'scaled:scale') is a p-vector of column
  # standard deviations of X if scale=TRUE, a p vector of ones
  # otherwise;
  attr(X, "scaled:center") = out$scaled_center
  attr(X, "scaled:scale")  = out$scaled_scale
  # 'attr(X,'d') is a p-vector of column sums of
  # X.standardized^2,' where X.standardized is the matrix X centered
  # by attr(X,'scaled:center') and scaled by attr(X,'scaled:scale').
  # Requires the package `matrixStats`.
  ## attr(X, "d") is applied when computing ERSS (`get_ER2()`), which is
  ## to estimate the residual variance.
  ## It seems that the squared X should not be weighted?!
  attr(X, "d") = out$d

  if (intercept) {  ## do not affect the estimation of beta (MLE)
    mean_y = mean(y)
    attr(y, "scaled:center") <-
      sapply(1:ncol(W), function(j) weighted.mean(x=y, w=W[,j]))
  }


  n <- nrow(X)
  p <- ncol(X)
  if (p > 1000 & !requireNamespace("Rfast",quietly = TRUE))
    warning_message("For an X with many columns, please consider installing",
                    "the Rfast package for more efficient credible set (CS)",
                    "calculations.", style='hint')


  # Check weight matrix W
  if (!is.null(W)){

    # Check whether the dimensions of W and X match
    if (!( (length(W) == nrow(X)) | all(dim(W) == dim(X)) ))
      stop("The dimensions of W and input X do not match!")
  }

  # Initialize clamp fit.
  s <- init_setup(n=n, p=p, maxL=maxL, family="linear",
                  scaled_prior_variance=scaled_prior_variance,
                  residual_variance=residual_variance,
                  prior_inclusion_prob=prior_inclusion_prob,
                  varY=as.numeric(var(y)),
                  standardize=standardize)
  s <- init_finalize(s)


  # Initialize elbo to NA.
  elbo = rep(as.numeric(NA),max_iter + 1)
  elbo[1] = -Inf;

  # Initialize log-likelihood into NA
  loglik = rep(as.numeric(NA),max_iter+1)  ## not required to know...
  loglik[1] = -Inf

  # (Expected) weighted sum of squared residual
  WRSS = rep(as.numeric(NA),max_iter+1)  ## not required to know...
  WRSS[1] = -Inf

  tracking = list()

  for (tt in 1:max_iter) {

    if (track_fit)
      tracking[[tt]] = clamp_slim(s)

    # if !is.null(seed), update the random seed in every iteration.
    if (!is.null(seed)) {seed <- seed + tt}

    s <- update_each_effect(X = X, y = y, s = s, W = W, model = model,
                          mle_estimator = mle_estimator,
                          mle_variance_estimator = mle_variance_estimator,
                          nboots = nboots,
                          seed = seed,
                          estimate_prior_variance = estimate_prior_variance,
                          estimate_prior_method = estimate_prior_method,
                          check_null_threshold = check_null_threshold,
                          abnormal_proportion = abnormal_proportion,
                          robust_method = robust_method,
                          robust_estimator = robust_estimator)

    # Compute objective before updating residual variance because part
    # of the objective s$kl has already been computed under the
    # residual variance before the update.
    # elbo[tt+1] = get_elbo(X, y, s, W)
    elbo[tt+1]   = get_elbo(s)
    WRSS[tt+1]   = mean(s$EWR2)
    loglik[tt+1] = mean(s$Eloglik)

    if (verbose){
      print(paste("#iteration:", tt,
                  "; objective:", round(elbo[tt+1], 2),
                  "; expected WRSS:", round(WRSS[tt+1], 2),
                  "; loglik:", round(loglik[tt+1], 2)))
    }

    ## convergence condition
    ## some outputs:
    if (mle_variance_estimator == "bootstrap") {

      ## If the approach selects nothing,
      ## then the objective would be a constant after a few iterations.
      if (abs(elbo[tt+1] - elbo[tt]) < tol) {
        s$converged = TRUE
        break
      }

      ## Otherwise, due to the randomness in bootstrapping,
      ## the ELBO may not converge to a specific value; instead, it may
      ## fluctuate around a certain range. Thus, we use the following
      ## after a "burn-in" period, if the average of `converge_window`
      ## successive ELBO converges, then break the loop.
        if (tt > (burn_in + converge_window + 1)) {
          if ( abs(mean(elbo[(tt - converge_window - 1) : (tt-1)]) -
                   mean(elbo[(tt - converge_window) : (tt)]) ) < tol ) {
              s$converged = TRUE
              break
          }
        }
    } else {
      if (abs(elbo[tt+1] - elbo[tt]) < tol) {
        s$converged = TRUE
        break
      }
    }





    if (estimate_residual_variance) {
      s$sigma2 <- pmax(residual_variance_lowerbound,
                       estimate_residual_variance_fun(X,y,s))

      if (s$sigma2 > residual_variance_upperbound)
        s$sigma2 <- residual_variance_upperbound
    }


  }

  # Remove first (infinite) entry, and trailing NAs.
  elbo = elbo[2:(tt+1)]
  s$elbo = elbo
  loglik = loglik[2:(tt+1)]
  s$loglik = loglik
  # Remove redundancies
  s$Eloglik = NULL
  s$EWR2 = NULL

  s$niter = tt

  if (is.null(s$converged)) {
    warning(paste("IBSS algorithm did not converge in",max_iter,"iterations!"))
    s$converged = FALSE
  }


  if (intercept) {
    # Estimate unshrunk intercept.
    # Treat this as an ordinary linear regression model
    # Here, the intercept and the estimated residual variance implicitly
    # represents the effect of the confounding factor on y.
    s$intercept <- mean_y -
      sum( colMeans(X) * (colSums(s$alpha * s$mu) / colSds(X)) )
    s$fitted <- s$Xr + mean_y
    # s$intercept <- mean_y -
    #   sum( attr(X,"scaled:center") *
    #          (colSums(s$alpha * s$mu)/attr(X,"scaled:scale")) )
  } else {
    s$intercept <- 0
    s$fitted <- s$Xr
  }


  s$fitted = drop(s$fitted)
  names(s$fitted) = `if`(is.null(names(y)),rownames(X),names(y))

  if (track_fit)
    s$trace = tracking

  # Credible Sets and PIPs
  if (!is.null(coverage) && !is.null(min_abs_corr)) {
    s$sets = clamp_get_cs(s,coverage = coverage,X = X,
                                  min_abs_corr = min_abs_corr,
                                  n_purity = n_purity)
    s$pip = clamp_get_pip(s,prune_by_cs = FALSE,prior_tol = prior_tol)
  }

  # (Re)name the outputs
  variable_names             <- colnames(X)
  names(s$pip)               <- variable_names
  colnames(s$alpha)          <- variable_names
  colnames(s$mu)             <- variable_names
  colnames(s$mu2)            <- variable_names
  colnames(s$betahat)        <- variable_names
  colnames(s$logBF_variable) <- variable_names

  # For prediction.
  # s$X_column_scale_factors  <- attr(X,"scaled:scale")
  # s$X_column_center_factors <- attr(X,"scaled:center")

  return(s)
}
