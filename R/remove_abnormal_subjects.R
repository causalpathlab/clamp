#' @title Remove abnormal points from the current matrix/array
#'
#' @param abnormal_index an array of indices of abnormal points.
#'
#' @param object a matrix or a vector that may contain abnormal points.
#' If \code{object} is a matrix, then \code{abnormal_index} indicates
#' the rows to be removed.
#'
#' @keywords internal
remove_abnormal_subjects <- function(abnormal_index, object) {
  if (!(is.matrix(object) | is.vector(object)))
    stop("Input object should be a matrix or a vector")

  abnormal_index <- unique(abnormal_index)

  if (is.matrix(object)) {
    len <- nrow(object)
  } else { # is.vector(object)
    len <- length(object)
  }

  if (!is.null(abnormal_index)) {
    subset_index <- (!(1:len) %in% abnormal_index)
  } else {
    subset_index <- rep(TRUE, times = len)
  }

  if (is.matrix(object)) {
    return(object[subset_index, , drop = F])
  } else {  # is.vector(object)
    return(object[subset_index])
  }

}
