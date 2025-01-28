# Compute the column means of X, the column standard deviations of X,
# and rowSums(Y^2), where Y is the centered and/or scaled version of
# X.
# W: an n by p weight matrix of the same size as X.
#
#' @importFrom Matrix rowSums
#' @importFrom Matrix colMeans
#' @importFrom matrixStats weightedSd
compute_colstats = function (X, W = NULL, center = TRUE, scale = TRUE) {
  n = nrow(X)
  p = ncol(X)

  if (is.null(W)) {
    # X is an ordinary dense or sparse matrix. Set sd = 1 when the
    # column has variance 0.
    if (center)
      scaled_center = colMeans(X,na.rm = TRUE)  ## `cm` in `susieR`
    else
      scaled_center = rep(0,p)
    if (scale) {
      scaled_scale = compute_colSds(X)
      scaled_scale[scaled_scale == 0] = 1
    } else
      scaled_scale = rep(1,p)

    # These two lines of code should give the same result as
    #
    #   Y = (t(X) - scaled_center)/scaled_scale
    #   d = rowSums(Y^2)

    # for all four combinations of "center" and "scale", but do so
    # without having to modify X, or create copies of X in memory. In
    # particular the first line should be equivalent to colSums(X^2).
    d = n * colMeans(X)^2 + (n-1) * compute_colSds(X)^2
    d = (d - n * scaled_center^2) / scaled_scale^2

  } else { # if (!is.null(W))

    # check the dimension of W
    if (!all(dim(W) == dim(X))) stop("The dimensions of W and X do not match!")

    # compute the weighted means and standard deviation.
    if (center)
      scaled_center = sapply(1:p, function(j) weighted.mean(x=X[,j], w=W[,j]))
    else
      scaled_center = rep(0, p)
    if (scale) {
      scaled_scale = sapply(1:p, function(j) weightedSd(x=X[,j], w=W[,j]))
      scaled_scale[scaled_scale == 0] = 1
    } else
      scaled_scale = rep(1, p)

    # column-wise standardized X^2
    # ...
    X_ = sweep(X, 2, scaled_center, "-")
    X_ = sweep(X_, 2, scaled_scale, "/")
    d = colSums(W * X_^2)

  }


  return(list(scaled_center = scaled_center,
              scaled_scale = scaled_scale,
              d = d))
}

# @title computes column standard deviations for any type of matrix
# @details This should give the same result as matrixStats::colSds(X),
#   but allows for sparse matrices as well as dense ones.
# @param X an n by p matrix of any type, e.g. sparse, dense.
# @return a p vector of column standard deviations.
#
#' @importFrom matrixStats colSds
#' @importFrom Matrix summary

compute_colSds = function(X) {
  if (is.matrix(X))
    return(colSds(X))
  else {  ## If X is not a ordinary matrix, i.e., X is a sparse matrix.
    n = nrow(X)
    X2 = apply_nonzeros(X,function (u) u^2)
    d = colMeans(X2) - colMeans(X)^2
    return(sqrt(d*n/(n-1)))
  }
}

# Apply operation f to all nonzeros of a sparse matrix.
#
#' @importFrom Matrix sparseMatrix
#'
apply_nonzeros <- function (X, f) {
  d <- summary(X)
  return(sparseMatrix(i = d$i,j = d$j, x = f(d$x), dims = dim(X)))
}

