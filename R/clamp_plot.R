#' @rdname clamp_plots
#'
#' @title Clamp Plots.
#'
#' @description \code{clamp_plot} produces a per-variable summary of
#'   the Clamp credible sets.
#'
#' @param model A SuSiE or Clamp fit, typically an output from
#'   \code{\link{clamp}} or one of its variants.
#'   For \code{gsuse_plot},
#'   the susie fit must have \code{model$z}(?), \code{model$PIP}, and may
#'   include \code{model$sets}. \code{model} may also be a vector of
#'   z-scores or PIPs.
#'
#' @param y A string indicating what to plot: either \code{"z_original"} for
#'   z-scores, \code{"z"} for z-score derived p-values on (base-10) log-scale,
#'   \code{"PIP"} for posterior inclusion probabilities,
#'   \code{"log10PIP"} for posterior inclusion probabiliities on the
#'   (base-10) log-scale. For any other setting, the data are plotted as
#'   is.
#'
#' @param include_intercept Boolean. If \code{include_intercept = FALSE},
#'   the intercept term is removed from the plot; the index of the intercept
#'   is specified in \code{intercept_index}.
#'
#' @param intercept_index a numeric number specifying the index of the intercept.
#'    It is required when \code{include_intercept=FALSE} (by default).
#'
#' @param effect_indices a numeric vector of indices of effect variables.
#'    Those points are highlighted in red.
#'
#' @param add_bar If \code{add_bar = TRUE}, add horizontal bar to
#'   signals in credible interval.
#'
#' @param pos This can be either be (1) a numeric vector of indices of
#'   subset of variables to plot, or (2) a list with the following list
#'   elements: \code{pos$attr}, \code{pos$start} and \code{pos$end},
#'   where \code{pos$attr} is a character string of the name of index
#'   variable in \code{model} object, and \code{pos$start} and
#'   \code{pos$end} are boundaries of indices to plot. See the provided
#'   examples.
#'
#' @param max_cs The largest credible set to display, either based on
#'   purity (set \code{max_cs} between 0 and 1), or based on size (set
#'   \code{max_cs > 1}).
#'
#' @param add_legend If \code{add_legend = TRUE}, add a legend to
#'   annotate the size and purity of each CS discovered. It can also be
#'   specified as location where legends should be added, e.g.,
#'   \code{add_legend = "bottomright"} (default location is
#'   \code{"topright"}).
#'
#' @param \dots Additional arguments passed to
#'   \code{\link[graphics]{plot}}.
#'
#' @return Invisibly returns \code{NULL}.
#'
#' @importFrom utils head
#' @importFrom stats pnorm
#' @importFrom graphics plot
#' @importFrom graphics segments
#' @importFrom graphics points
#' @importFrom graphics legend
#' @importFrom graphics par
#'
#' @export
#'
clamp_plot <- function (model, y,
                         include_intercept = FALSE,
                         intercept_index = NULL,
                         effect_indices = NULL,
                         add_bar = FALSE,
                         pos = NULL,
                         max_cs = 400, add_legend = NULL, ...) {
  is_clamp <- inherits(model,c("susie", "clamp"))
  ylab <- y
  color <- c(
    "dodgerblue2",
    "green4",
    "#6A3D9A", # purple
    "#FF7F00", # orange
    "gold1",
    "skyblue2", "#FB9A99", # lt pink
    "palegreen2",
    "#CAB2D6", # lt purple
    "#FDBF6F", # lt orange
    "gray70", "khaki2",
    "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
    "darkturquoise", "green1", "yellow4", "yellow3",
    "darkorange4", "brown"
  )

  if (y == "PIP"){
    if (is_clamp)
      yy <- model$pip
    else
      yy <- model
  } else if (y == "log10PIP") {
    if (is_clamp)
      yy <- log10(model$pip)
    else
      yy <- log10(model)
    ylab <- "log10(PIP)"
  } else {
    if (is_clamp){
      stop("Please specify y for the plot.")
    }
    yy <- model
  }

  if (!include_intercept) {
    if (class(model) == "clamp") {
      if (!is.null(model$family) & (model$family %in% c("logistic", "poisson"))) {
        intercept_index <- length(model$pip)
        yy <- yy[-intercept_index]
      }
    }
  }

  if(is.null(pos))
    pos <- 1:length(yy)
  start <- 0
  if (inherits(pos, "list")) {

    # Check input.
    if (is.null(pos$attr) || is.null(pos$start) || is.null(pos$end))
      stop("pos argument should be a list of list(attr=,start=,end=)")
    if (!(pos$attr %in% names(model)))
      stop(paste("Cannot find attribute",pos$attr,"in input model object"))
    if (pos$start >= pos$end)
      stop("Position start should be smaller than end")
    start <- min(min(model[[pos$attr]]),pos$start)
    end <- max(max(model[[pos$attr]]),pos$end)

    # Add zeros to alpha and p.
    alpha <- matrix(0,nrow(model$alpha),end - start + 1)
    new_yy <- rep(min(yy),end - start + 1)
    pos_with_value <- model[[pos$attr]] - start + 1
    new_yy[pos_with_value] <- yy
    alpha[,pos_with_value] <- model$alpha
    yy <- new_yy
    model$alpha <- alpha

    # Adjust model$cs.
    if (!is.null(model$sets$cs)) {
      for (i in 1:length(model$sets$cs))
        model$sets$cs[[i]] <- pos_with_value[model$sets$cs[[i]]]
    }

    # Change "pos" object to be indices.
    start_adj <- -min(min(model[[pos$attr]]) - pos$start,0)
    end_adj <- max(max(model[[pos$attr]]) - pos$end,0)
    pos <- (1 + start_adj):(length(yy) - end_adj)
  } else {
    if (!all(pos %in% 1:length(yy)))
      stop("Provided position is outside the range of variables")
  }
  legend_text <- list(col = vector(), purity = vector(), size = vector())
  # scipen0 = options()$scipen
  # options(scipen = 10)
  args <- list(...)
  if (!exists("xlab", args)) args$xlab <- "variable"
  if (!exists("ylab", args)) args$ylab <- ylab
  if (!exists("pch", args)) args$pch <- 16
  args$x <- pos + start
  args$y <- yy[pos]
  do.call(plot, args)
  if (is_clamp && !is.null(model$sets$cs)) {
    for(i in rev(1:nrow(model$alpha))){
      if (!is.null(model$sets$cs_index) && !(i %in% model$sets$cs_index))
        next
      purity <- model$sets$purity[which(model$sets$cs_index == i),1]
      if (!is.null(model$sets$purity) && max_cs < 1 && purity >= max_cs) {
        x0 <- intersect(pos,model$sets$cs[[which(model$sets$cs_index == i)]])
        y1 <- yy[x0]
      } else if (n_in_CS(model, model$sets$requested_coverage)[i] < max_cs) {
        x0 <- intersect(pos,
                        which(in_CS(model,model$sets$requested_coverage)[i,] > 0))
        y1 <- yy[x0]
      } else {
        x0 <- NULL
        y1 <- NULL
      }
      if (is.null(x0))
        next
      if (add_bar) {
        y0 <- rep(0,length(x0))
        x1 <- x0
        segments(x0+start,y0,x1+start,y1,lwd = 1.5,col = "gray")
      }
      points(x0+start,y1,col = head(color,1),cex = 1.5,lwd = 2.5)
      legend_text$col <- append(head(color,1), legend_text$col)

      # Rotate color.
      color <- c(color[-1],color[1])
      legend_text$purity <- append(round(purity,4),legend_text$purity)
      legend_text$size <- append(length(x0),legend_text$size)
    }
    if (length(legend_text$col) > 0 && !is.null(add_legend) &&
        !identical(add_legend, FALSE)) {

      # Plot legend.
      text <- vector()
      for (i in 1:length(legend_text$col)) {
        if (legend_text$size[i] == 1)
          text[i] <- paste0("L",i,": C=1")
        else
          text[i] <- paste0("L",i,": C=",legend_text$size[i],"/R=",
                            legend_text$purity[i])
      }
      if (!(add_legend %in% c("bottomright", "bottom", "bottomleft", "left",
                              "topleft", "top", "topright", "right", "center"))) {
        add_legend <- "topright"
      }
      legend(add_legend,text,bty = "n",col = legend_text$col,cex = 0.65,
             pch = 15)
    }
  }
  if (!is.null(effect_indices))
    points(pos[effect_indices] + start, yy[effect_indices] + start,
           col = 2,pch = 16)
  # options(scipen = scipen0)
  return(invisible())
}
