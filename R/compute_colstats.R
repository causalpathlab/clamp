# Compute the column means of X, the column standard deviations of X,
# and rowSums(Y^2), where Y is the centered and/or scaled version of
# X.
#
#' @importFrom Matrix rowSums
#' @importFrom Matrix colMeans
compute_colstats = function (X, center = TRUE, scale = TRUE) {
  n = nrow(X)
  p = ncol(X)

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
  #   Y = (t(X) - cm)/csd
  #   d = rowSums(Y^2)

  # for all four combinations of "center" and "scale", but do so
  # without having to modify X, or create copies of X in memory. In
  # particular the first line should be equivalent to colSums(X^2).
  d = n * colMeans(X)^2 + (n-1) * compute_colSds(X)^2
  d = (d - n * scaled_center^2) / scaled_scale^2

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

