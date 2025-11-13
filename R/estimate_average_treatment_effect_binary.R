#' @rdname estimate_average_treatment_effect
#'
#' @title Estimate the average treatment effect (ATE) via inverse probability weighting (IPW)
#'
#' @description
#' This function provides two approaches to provides IPW estimates of the ATE.
#' One IPW estimator is the modified Horvitz-Thompson (mHT) estimator, and
#' the alternative is the weighted least-squares (WLS) estimator. When the input
#' \code{\mathbf{X}} is binary, the IPW estimates from these two methods should
#' be equal.
#'
#' @param y An (n by 1) vector of responses.
#'
#' @param X An (n by p) data matrix; each entry would be 0 or 1.
#'
#' @param W An (n by p) weight matrix.
#'
#' @param mle_estimator Estimation approach. \code{"mHT"} stands for the
#' modified Horvitz-Thompson estimator. \code{"WLS"} stands for the weighted
#' least-squares estimator.
#'
#' @importFrom matrixStats weightedSd
#' @importFrom stats weighted.mean
#'
#' @keywords internal
estimate_average_treatment_effect_binary <- function(X, y, W,
                                              mle_estimator = c("mHT", "WLS"),
                                              standardize = NULL,
                                              centralize = NULL) {


  if (!all(dim(X) == dim(W))) stop("Dimensions of X an W do not match!")
  if (nrow(X) != length(y))   stop("length(y) != nrow(X)")

  mle_estimator <- match.arg(mle_estimator)

  switch(mle_estimator,

         "mHT" = {

           # When applying the modified Horvitz-Thompson estimator
           # the input X and response Y do not need to be scaled.
           # Each entry of X should be either 0 or 1.
           deltahat <-  colSums( sweep(W*X, 1, y, "*") ) / colSums( W*X ) -
             colSums( sweep(W*(1-X), 1, y, "*") ) / colSums( W*(1-X) )
         },

         "WLS" = {

           # Scale X:
           ## 1. When inputting the original data matrix `X`,
           ## `attr(X, "scaled:center)` and `attr(X, "scaled:scale")` would
           ## not be NULL. In this case, we centralize and standardize each
           ## column of X according to their coresponding center and scale.
           ## 2. If inputing the bootstrap data matrix `Xboot` as `X`,
           ## `attr(X, "scaled:center)` and `attr(X, "scaled:scale")` would
           ## be NULL. In this case, whether each column of `Xboot` are
           ## centralized and standardized depends on the arguments
           ## `centralize` and `standardized`.

           if ( !is.null(attr(X, "scaled:center")) &
                !is.null(attr(X, "scaled:scale"))) {
             X_ <- sweep(X, 2, attr(X, "scaled:center"), "-")
             X_ <- sweep(X_, 2, attr(X, "scaled:scale"), "/")
           }
           else {  ## If bootstrapping

             if (centralize) {
               X_ <- sapply(1:ncol(W),
                            function(j) X[,j] - weighted.mean(x=X[,j], w=W[,j]))
             }
             else {
               X_ <- X
             }

             if (standardize) {
               X_ <- sapply(1:ncol(W),
                            function(j) X_[,j] / weightedSd(x=X[,j], w=W[,j]))
             }
             else {
               X_ <- X_
             }

           }

           # Scale y: y_ is a n by p matrix.
           ## 1. When inputting the original response vector `y`,
           ## `attr(y, "scaled:center")` is not NULL. In this case, we create
           ## an n by p matrix `y_`, of which each column is a centralized `y`
           ## with respect to the weighted mean `weighted.mean(y, W[,j])`.
           ## 2. When inputting the bootstrap response vector `yboot` as `y`,
           ## `attr(y, "scaled:center")`
           if ( !is.null(attr(y, "scaled:center")) ) {
             y_ <- outer(c(y), attr(y, "scaled:center"), "-")
           } else {
             if (centralize) {
               y_ <- sapply(1:ncol(W), function(j) y - weighted.mean(y, W[,j]))
             } else {
               y_ <- matrix(rep(y, times = ncol(X_)), ncol = ncol(X_))
             }
           }

           # Obtain weighted least-squares estimates
           wxy <- colSums(W * X_ * y_)
           wx2  <- colSums(W * X_^2)
           deltahat <- wxy / wx2
         }
         )

  # Return
  if ( !is.null(attr(X, "scaled:scale")) )  ## if not bootstrapping
    return( deltahat/attr(X, "scaled:scale") )
  else   ## if bootstrapping, return the estimates directly.
    return(deltahat)

}
