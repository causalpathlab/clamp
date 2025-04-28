#' @rdname estimate_average_treatment_effect_categorical
#'
#' @title Estimate the average treatment effect via inverse probability weighting (IPW)
#'
#' @description
#' This function computes the IPW estimates of the average treatment effect.
#' One IPW estimator is the modified Horvitz-Thompson (mHT) estimator.
#'
#' @param y An (n by 1) vector of responses.
#'
#' @param X An (n by p) data matrix; each entry would be 0 or 1. Each column of
#' X should be named in the format of \code{variable_level}.
#'
#' @param W An (n by p) weight matrix.
#'
#' @param causal_effect_estimator Estimation approach. \code{"mHT"} stands for
#' the modified Horvitz-Thompson estimator.
#'
#' @returns a list of the following items:
#'
#' \item{thetahat0}{a named vector of average effect estimates of
#'  baseline levels.}
#'  \item{deltahat}{a named vector of average treatment effect estimates.}
#'
#' @importFrom matrixStats weightedSd  ## not used
#' @importFrom stats weighted.mean     ## not used
#' @importFrom reshape2 dcast
#' @importFrom reshape2 melt
#'
#' @keywords internal
estimate_average_treatment_effect_categorical <- function(X, y, W,
                                           causal_effect_estimator = c("mHT")
                                           ) {


  if (!all(dim(X) == dim(W))) stop("Dimensions of X an W do not match!")
  if (nrow(X) != length(y))   stop("length(y) != nrow(X)")

  causal_effect_estimator <- match.arg(causal_effect_estimator)

  if (is.null(colnames(X)))
    stop("Columns of X should be named in the format of `variable_level`.")

  if (causal_effect_estimator == "mHT") {
    # When applying the modified Horvitz-Thompson estimator
    # the input X and response Y do not need to be scaled.
    # Each entry of X should be either 0 or 1.

    # the modified Horvitz-Thompson estimator for average causal effect
    # 1 by ncol(X) named matrix
    thetahat <- colSums( sweep(W*X, 1, y, "*") ) / colSums( W*X )

    if (is.null(names(thetahat)))
      stop("Elements of thetahat should be named in the format of `variable_level`.")

    ## Compute deltahat
    ### 1. Extract variable and level from names
    var_level <- strsplit(names(thetahat), "_")
    vars <- sapply(var_level, `[`, 1)
    levels <- sapply(var_level, `[`, 2)

    ### 2. Get level 0 (baseline level) values for each variable
    thetahat0 <- thetahat[levels == "0"]
    names(thetahat0) <- vars[levels == "0"]

    ### 3. Compute differences for non-baseline levels
    non_baseline <- levels != "0"
    deltahat <- thetahat[non_baseline] - thetahat0[vars[non_baseline]]
    names(deltahat) <- names(thetahat)[non_baseline]
  }

  return(list(
    deltahat = deltahat,
    thetahat0 = thetahat0
  ))
}



# estimate the average marginal potential outcomes (po) using IPW methods
#' @importFrom reshape2 dcast
#' @importFrom reshape2 melt
#'
#' @param y An (n by 1) vector of responses.
#'
#' @param X An (n by p?) data matrix; each entry would be 0 or 1. Each column of
#' X should be named in the format of \code{variable_level}.
#'
#' @param W An (n by p?) weight matrix.
#'
#' @param thetahat0 a vector of estimated baseline effects
#'
#' @param deltahat a vector of estimated treatment effects.
#'
#' @returns a vector of estimated average potential outcomes corresponding to
#' each margin, in which each entry is named in the format of "variable_level".
#'
#' @keywords internal
estimate_average_marginal_po <- function(X = NULL, W = NULL, y = NULL,
                                         thetahat0 = NULL, deltahat = NULL,
                                         .sort = TRUE,
                                         .variable_prefix = "X") {

  if ( !(is.null(X) | is.null(W) | is.null(y)) ) {
    thetahat <- colSums( sweep(W*X, 1, y, "*") ) / colSums( W*X )

  } else {
    ## retrieve the estimated causal effects from thetahat0 and deltahat.
    if (is.null(thetahat0) | is.null(deltahat)) {
      stop("Please specify both thetahat0 and deltahat.")
    }

    ### 1. Get variable names back from deltahat
    var_level <- strsplit(names(deltahat), "_")
    vars <- sapply(var_level, "[", 1)

    ### 2. Reconstruct the non-baseline levels
    non_baseline_values <- deltahat + thetahat0[vars]
    names(non_baseline_values) <- names(deltahat)

    ### 3. Add back the baseline (level 0) entries
    thetahat <- c(thetahat0, non_baseline_values)
    names(thetahat)[1:length(thetahat0)] <- paste0(names(thetahat0), "_0")

    ### 4. Sort thetahat based on their variable names (optional)
    if (.sort) {
      # Extract variable and level from names
      var_level <- do.call(rbind, strsplit(names(thetahat), "_"))
      vars <- var_level[,1]
      levels <- as.numeric(var_level[,2])

      # Order by variable, then level
      var_num <- as.numeric(sub(.variable_prefix, "", vars))
      ## so the variable name should begin with `.variable_prefix`

      # Get the ordering index
      ord <- order(var_num, levels)

      # Reorder the vector
      thetahat <- thetahat[ord]
    }
  }

  return(thetahat)
}


# if (TRUE) {
#   # average treatment effect estimator, deltahat
#
#   # suppose there is a vector indicating the columns of X.
#   # suppose the column names are separated by "_".
#   col_indices_df <- as.data.frame(
#     do.call(rbind, strsplit(names(thetahat), "_")))
#   # 1st column: variable; 2nd column: level
#   colnames(col_indices_df) <- c("variable", "level")
#
#   thetahat_df <- cbind(col_indices_df, thetahat = thetahat)
#   thetahat_df <- dcast(thetahat_df, level ~ variable,
#                        value.var = "thetahat")
#   rownames(thetahat_df) <- thetahat_df[,"level"]
#   thetahat_df <- thetahat_df[, colnames(thetahat_df) != "level"]
#   deltahat_df <-
#     sweep(as.matrix(thetahat_df[-1, colnames(thetahat_df) != "level"]), 2,
#           as.matrix(thetahat_df[1,  colnames(thetahat_df) != "level"]), "-")
#   # each column represents a variable, and each row represents a level
#   deltahat_df <- melt(deltahat_df, varnames = c("level", "variable"),
#                       value.name = "deltahat")
#   deltahat <- deltahat_df[,"deltahat"]
#   names(deltahat) <- paste0(deltahat_df[,"variable"], "_",
#                             deltahat_df[,"level"])
#   # remove empty levels (NA)
#   deltahat <- deltahat[!is.na(deltahat)]
# }

# coefs <- append(thetahat0, deltahat)
# col_indices_df <- as.data.frame(
#   do.call(rbind, strsplit(names(coefs), "_")))
# # 1st column: variable; 2nd column: level
# colnames(col_indices_df) <- c("variable", "level")
#
# coefs <- cbind(col_indices_df, coefs)
# coefs <- dcast(coefs, level ~ variable, value.var = "coefs")
# coefs[-1, colnames(coefs)!= "level"] <-
#   sweep(as.matrix(coefs[-1, colnames(coefs)!= "level"]), 2,
#         as.matrix(coefs[1,  colnames(coefs)!= "level"]), "+")
# thetahat_df <- melt(coefs, varnames = c("level", "variable"),
#                     value.name = "thetahat")
# thetahat <- thetahat_df[,"thetahat"]
# names(thetahat) <- paste0(thetahat_df[, "variable"], "_",
#                           thetahat_df[, "level"])
# # remove empty levels (NA)
# thetahat <- thetahat[!is.na(thetahat)]

