#' @title Summarize posterior estimations of regression coefficients
#'
#' @description
#' This function returns, by default, \code{top_n=10} coefficients with the
#' highest PIPs, sorted in decreasing order.
#'
#' @param object a clamp fit
#'
#' @param subset_variables an array of numeric indices or variable names
#'    presented in the output. If \code{level_names = NULL},
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
#' @param credible_interval logical, whether to return equal-tailed
#'    credible intervals from the normal variational distributions.
#'
#' @param coverage numeric, between 0 and 1, the coverage probability
#'    of credible intervals.
#'
#' @param digits integer indicating the number of decimal places to be used.
#'    if \code{digits == NULL}, the outputput will not be rounded.
#'
#' @returns The function outputputs a data.frame with \code{top_n} variables with
#'    the highest PIPs, or variables specified in \code{variables}.
#'    Variables are sorted in decreasing order of PIP values.
#'    The data.frame contains \emph{variable name}, \emph{PIP},
#'    \emph{posterior mean} (\code{posterior_mean}),
#'    \emph{posterior sd} (\code{posterior_sd}),
#'    and, if \code{cred_int=TRUE}, \emph{credible interval}
#'    (\code{crI_lower} and \code{crI_upper})
#'    of each variable.
#'
#' @rdname clamp_summarize_coefficients
#'
#' @importFrom stats qnorm
#'
#' @export
#'
clamp_summarize_coefficients <- function(object,
                                         subset_variables = NULL,
                                         decreasing = TRUE,
                                         top_n = 10,
                                         credible_interval = TRUE,
                                         coverage = 0.95,
                                         digits = 4) {

  if (is.null(object$level_pip) || is.null(object$mu) || is.null(object$mu2))
    stop("Either PIP, mu, or mu2 is not available.")

  # extract variable names
  if (is.null(names(object$level_pip))) {
    level_names <- 1:length(object$level_pip)
    variable_names <- 1:length(object$variable_pip)
  } else {
    level_names <- names(object$level_pip)
    variable_names <- sapply(strsplit(names(object$level_pip), "_"), `[`, 1)
  }

  output <- data.frame(
    variable  = variable_names,
    level     = level_names,
    level_pip = object$level_pip,
    posterior_mean  = clamp_get_posterior_mean(object),
    posterior_sd    = clamp_get_posterior_sd(object),
    variable_pip = object$variable_pip[variable_names]
  )

  if (credible_interval) {
    if (!is.numeric(coverage) | coverage < 0 | coverage > 1) {
      stop("Input probability should between 0 and 1")
    }

    output$crI_lower <- output$posterior_mean +
      qnorm((1 - coverage)/2) * output$posterior_sd
    output$crI_upper <- output$posterior_mean -
      qnorm((1 - coverage)/2) * output$posterior_sd
  }

  output <- output[order( output$level_pip,
                          output$variable_pip,
                          abs(output$posterior_mean),
                          decreasing = decreasing ), , drop=F]

  ## round the numeric in the output
  if (!is.null(digits)) {
    is_numeric_cols <- sapply(output, is.numeric)
    output[, is_numeric_cols] <-
      round(output[, is_numeric_cols, drop = F], digits = digits)
  }

  if (is.null(subset_variables)) {
    num_variables_output <- min(top_n, length(object$variable_pip))
    num_levels_output <- num_variables_output *
      ( length(object$level_pip) %/% length(object$variable_pip) )

    print.data.frame(
      output[1:num_levels_output, , drop=F], row.names = F
      )

  } else {
    print.data.frame(
      output[output[,"variable"] %in% subset_variables, , drop=F],
      row.names = F
      )
  }

  # return(output)
}
