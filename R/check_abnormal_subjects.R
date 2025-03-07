#' @title Sensitivity check
#'
#' @description
#' In iterative weighted least square methods,
#' each data point is transformed into a corresponding pseudo-response with
#' a corresponding pseudo-variance.
#' These transformations would inevitably point to some abnormal subjects
#' whose pseudo-response or pseudo-variance is an abnormal value,
#' such as Inf/-Inf and NaN(?).
#' We take these subjects with abnormal pseudo-response/pseudo-variance
#' as \dQuote{abnormal subjects}.
#' \code{check_abnormal_subjects()} return the indices of abnormal
#' subjects (if any).
#'
#' @param values A vector of length n, to check whether if it contains any
#' \code{Inf} or \code{NAN}.
#'
#' @returns index vector named \code{abnormal_index} such that
#'  \code{values[abnormal_index]} is either \code{Inf} or \code{NAN}.

check_abnormal_subjects <- function(values) {
  # check if any infinite
  abnormal_index <- which(is.infinite(values) | is.nan(values))
  if (length(abnormal_index) == 0) abnormal_index <- NULL

  return(abnormal_index = abnormal_index)
}
