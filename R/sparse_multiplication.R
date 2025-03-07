# @title Computes standardized.X %*% b using sparse multiplication trick
# Copied and modified from susieR.
#
# @param X an n by p unstandardized matrix with three attributes:
# attr(X,"scaled:center"), attr(X,"scaled:scale") and attr(X,"d")
# @param b a p-vector
# @return an n-vector
#
#' @importFrom Matrix t
#' @importFrom Matrix tcrossprod
compute_Xb = function (X, b) {

  # scaled_center = attr(X,"scaled:center")
  # scaled_scale = attr(X,"scaled:scale")

  # Scale Xb.

  # When X is an ordinary sparse/dense matrix.
  # scaled.Xb = tcrossprod(X,t(b/scaled_scale))  ## if X itself is standardized
  Xb = tcrossprod(X,t(b))  ## if X itself is not standardized.

  # Center Xb.
  # Xb = scaled.Xb - sum(scaled_center*b/scaled_scale)
  # This part is removed because in `clamp()` and `gsusie()`
  # the mean of y is not deduced, i.e., y = y - mean(y) is no longer operatred.

  return(as.numeric(Xb))
}


# # @title Computes t(standardized.X) %*% y using sparse multiplication trick
# # @param X an n by p unstandardized matrix with three attributes:
# # attr(X,"scaled:center"), attr(X,"scaled:scale") and attr(X,"d")
# # @param y an n vector
# # @return a p vector
# #
# #' @importFrom Matrix t
# #' @importFrom Matrix crossprod
# compute_Xty = function (X, y) {
#   scaled_center = attr(X,"scaled:center")
#   scaled_scale = attr(X,"scaled:scale")
#   ytX = crossprod(y,X)
#
#   # Scale Xty.
#
#   # When X is an ordinary sparse/dense matrix.
#   scaled.Xty = t(ytX/scaled_scale)
#
#   # Center Xty.
#   centered.scaled.Xty = scaled.Xty - scaled_center/scaled_scale * sum(y)
#   return(as.numeric(centered.scaled.Xty))
# }

# @title Computes M %* %t(standardized.X) using sparse multiplication trick
# @param M a L by p matrix
# @param X an n by p unstandardized matrix with three attributes:
# attr(X,"scaled:center"), attr(X,"scaled:scale") and attr(X,"d")
# @return a L by n matrix
#
#' @importFrom Matrix t
compute_MXt = function (M, X) {

  # scaled_center = attr(X, "scaled:center")
  # if (is.null(scaled_center)) scaled_center = rep(0, times = ncol(X))
  #default: no centralization

  # scaled_scale = attr(X, "scaled:scale")
  # if (is.null(scaled_scale)) scaled_scale = rep(1, times = ncol(X))
  #default: no standardization

  # When X is an ordinary sparse/dense matrix.
  ## De-centralizing part is removed because in `clamp()` and `gsusie()`
  ## the mean of y is not deduced, i.e., y = y - mean(y) is no longer operated.
  ## Besides, de-standardizing the coefficient part is also removed because
  ## the mean is also
  # return(as.matrix(t(X %*% (t(M)/scaled_scale)) -
  #                    drop(M %*% (scaled_center/scaled_scale))))
  return( as.matrix(t(X %*% (t(M)))) )
}

