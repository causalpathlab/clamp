#' @title Summarize posterior estimations of regression coefficients
#'
#' @description
#' This function returns, by default, \code{top_n=10} coefficients with the
#' highest PIPs, sorted in decreasing order.
#'
#' @param object a clamp fit
#'
#' @param subset_vars an array of numeric indices or variable names
#'    presented in the output. If \code{var_names = NULL},
#'    by default the output contains coefficients of all variables
#'    if \code{ncol(X)<= 10} or \code{top_n} variables with highest PIPs
#'    if \code{ncol(X) > 10}.
#'
#' @param decreasing logical, whether to list the variables
#'    in decreasing order of PIP values.
#'
#' @param top_n numeric, the number of coefficients to be displayed. By default
#' \code{top_n=min(10, length(object$pip))}.
#'
#' @param cred_int logical, whether to return equal-tailed credible intervals
#'
#' @param coverage numeric, between 0 and 1, the coverage probability
#'    of credible intervals.
#'
#' @param digits integer indicating the number of decimal places to be used.
#'    if \code{digits == NULL}, the output will not be rounded.
#'
#' @returns The function outputs a data.frame with \code{top_n} variables with
#'    the highest PIPs, or variables specified in \code{vars}.
#'    Variables are sorted in decreasing order of PIP values.
#'    The data.frame contains \emph{variable name}, \emph{PIP},
#'    \emph{posterior mean} (\code{PostMean}),
#'    \emph{posterior sd} (\code{PostSD}),
#'    and, if \code{cred_int=TRUE}, \emph{credible interval}
#'    (\code{CI_lower} and \code{CI_upper})
#'    of each variable.
#'
#' @rdname summarize_coefficients
#'
#' @importFrom stats qnorm
#'
#' @export
#'
summarize_coefficients <- function(object, subset_vars = NULL,
                                      decreasing = TRUE,
                                      top_n = 10,
                                      cred_int = TRUE,
                                      coverage = 0.95,
                                      digits = 4) {

  if (is.null(object$pip) || is.null(object$mu) || is.null(object$mu2))
    stop("Either PIP, mu, or mu2 is not available.")

  # extract variable names
  if (is.null(names(object$pip))) {
    var_names <- 1:length(object$pip)
  } else {
    var_names <- names(object$pip)
  }

  out <- data.frame(
    variable  = var_names,
    PIP       = object$pip,
    PostMean  = clamp_get_posterior_mean(object),
    PostSD    = clamp_get_posterior_sd(object)
  )

  if (cred_int) {
    if (!is.numeric(coverage) | coverage < 0 | coverage > 1) {
      stop("Input probability should between 0 and 1")
    }

    out$CrI_lower <- out$PostMean + qnorm((1 - coverage)/2) * out$PostSD
    out$CrI_upper <- out$PostMean - qnorm((1 - coverage)/2) * out$PostSD
  }

  # rownames(out) <- var_names

  out <- out[order(out$PIP, decreasing = decreasing), ]
  if (is.null(subset_vars)) {
    if (!is.null(digits)) {
      print.data.frame(cbind(variable = out[1:min(top_n, length(object$pip)), 1],
                             round(out[1:min(top_n, length(object$pip)),
                                       2:ncol(out)],
                                   digits)), row.names = F)
    } else {
      print.data.frame(out[1:min(top_n, length(object$pip)), ], row.names = F)
    }

  } else {
    if (!is.null(digits)) {
      print.data.frame(cbind(variable = out[out[, 1] %in% subset_vars, 1],
                             round(out[out[, 1] %in% subset_vars,
                                       2:ncol(out)],
                                   digits)), row.names = F)
    } else {
      print.data.frame(out[out[, 1] %in% subset_vars, ], row.names = F)
    }
  }

  # return(out)
}
