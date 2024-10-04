#' @title Summarize Susie Fit.
#' Copied from susieR.
#'
#' @description \code{summary} method for the \dQuote{clamp} class.
#'
#' @param object A clamp fit.
#'
#' @param \dots Additional arguments passed to the generic \code{summary}
#'   or \code{print.summary} method.
#'
#' @return \code{summary.clamp} returns a list containing a data frame
#'   of variables and a data frame of credible sets.
#'
#' @method summary clamp
#'
#' @export summary.clamp
#'
#' @export
#'
summary.clamp = function (object, ...) {
  if (is.null(object$sets))
    stop("Cannot summarize Clamp object because credible set information ",
         "is not available")
  variables = data.frame(cbind(1:length(object$pip),object$pip,-1))
  colnames(variables) = c("variable","variable_prob","cs")
  rownames(variables) = NULL
  if (!is.null(object$sets$cs)) {
    cs = data.frame(matrix(NA,length(object$sets$cs),5))
    colnames(cs) = c("cs","cs_logBF","cs_avg_r2","cs_min_r2","variable")
    for (i in 1:length(object$sets$cs)) {
      variables$cs[variables$variable %in% object$sets$cs[[i]]] =
        object$sets$cs_index[[i]]
      cs$cs[i] = object$sets$cs_index[[i]]
      cs$cs_logBF[i]  = object$logBF[cs$cs[i]]
      cs$cs_avg_r2[i] = object$sets$purity$mean.abs.corr[i]^2
      cs$cs_min_r2[i] = object$sets$purity$min.abs.corr[i]^2
      cs$variable[i] = paste(object$sets$cs[[i]],collapse=",")
    }
    variables = variables[order(variables$variable_prob,decreasing = TRUE),]
  } else
    cs = NULL
  out = list(vars = variables,cs = cs)
  class(out) = c("summary.clamp","list")
  return(out)
}

#' @rdname summary.clamp
#'
#' @param x A clamp summary.
#'
#' @method print summary.clamp
#'
#' @export print.summary.clamp
#'
#' @export
#'
print.summary.clamp = function (x, ...) {
  cat("\nVariables in credible sets:\n\n")
  print.data.frame(x$vars[which(x$vars$cs > 0),],row.names = FALSE)
  cat("\nCredible sets summary:\n\n")
  print.data.frame(x$cs,row.names = FALSE)
}
