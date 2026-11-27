#' Plot share ratio elasticities
#'
#' Draws the share ratio elasticities computed by [ShareRatioElasticities()] as
#' a matrix of points, with the share ratios of \eqn{Y} in columns and those of
#' \eqn{X} in rows, as in Figure 15 of Dargel and Thomas-Agnan (2024).
#'
#' @details
#' The elasticities are grouped into `n_bins` bins of width `width`. Each value
#' goes to the bin with the closest centre. A value outside the range covered by
#' the bins is put in the nearest outer bin.
#'
#' Elasticities that are not finite (when the share ratio of \eqn{X} does not
#' change in the direction `Xdir`) are not drawn.
#'
#' @param x an object of class "ShareRatioElasticities"
#' @param n_bins an odd integer of at least 3, the number of bins
#' @param width a positive number, the width of the bins. If `NULL` (the
#'   default), the width is `max(abs(x$SRE)) / (n_bins / 2)`, so that the
#'   largest absolute value falls in an outer bin
#' @param max_size a positive number, the size of the largest symbols in mm
#' @param ... not used
#' @return a ggplot object
#' @author
#'   - Lukas Dargel
#'   - Anh-Dung Le
#' @references
#'   - Dargel, Lukas and Christine Thomas-Agnan, “Pairwise share ratio interpretations of compositional regression models”, Computational Statistics & Data Analysis 195 (2024), p. 107945
#' @importFrom grDevices hcl.colors
#' @exportS3Method
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   res <- lmCoDa(
#'     ilr(cbind(left, right, extreme_right)) ~
#'       ilr(cbind(Age_1839, Age_4064, Age_65plus)),
#'     data = election
#'   )
#'   sre <- ShareRatioElasticities(
#'     res,
#'     Xvar = "cbind(Age_1839, Age_4064, Age_65plus)"
#'   )
#'   plot(sre, n_bins = 7)
#' }
plot.ShareRatioElasticities <- function(
  x,
  n_bins,
  width = NULL,
  max_size = 5,
  ...
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The package ggplot2 is needed to plot share ratio elasticities!")
  }
  if (!all(c("Yj", "Yl", "Xj", "Xl", "SRE") %in% names(x))) {
    stop("x must have the columns Yj, Yl, Xj, Xl and SRE!")
  }
  if (missing(n_bins)) {
    stop("n_bins must be given!")
  }
  if (
    !(is.numeric(n_bins) &&
      length(n_bins) == 1 &&
      !is.na(n_bins) &&
      n_bins >= 3 &&
      n_bins %% 2 == 1)
  ) {
    stop("n_bins must be an odd integer of at least 3!")
  }
  if (
    !(is.null(width) ||
      (is.numeric(width) && length(width) == 1 && isTRUE(width > 0)))
  ) {
    stop("width must be NULL or a positive number!")
  }
  if (
    !(is.numeric(max_size) && length(max_size) == 1 && isTRUE(max_size > 0))
  ) {
    stop("max_size must be a positive number!")
  }

  x <- as.data.frame(x)
  finite <- is.finite(x$SRE)
  if (!all(finite)) {
    message(
      sum(!finite),
      " share ratio elasticities are not finite and are not plotted."
    )
    x <- x[finite, ]
  }
  if (nrow(x) == 0) {
    stop("There is no finite share ratio elasticity to plot!")
  }

  # bins: centres k * width for k = -half, ..., half
  half <- (n_bins - 1) / 2
  if (is.null(width)) {
    width <- max(abs(x$SRE)) / (n_bins / 2)
    if (width == 0) {
      width <- 1
    }
  }
  centers <- seq(-half, half) * width
  bin_labels <- format(signif(centers, 3), trim = TRUE)

  # the closest centre of each value; a value on an edge goes to the lower bin,
  # values beyond the outer edges go to the outer bins
  closest_center <- function(v) {
    for (cen in centers) {
      if (v <= cen + width / 2) {
        return(cen)
      }
    }
    max(centers)
  }
  x$SREB <- factor(
    vapply(x$SRE, closest_center, numeric(1)),
    levels = centers,
    labels = bin_labels
  )

  # share ratios as axis labels, grouped by denominator
  y_parts <- unique(c(x$Yl, x$Yj))
  x_parts <- unique(c(x$Xl, x$Xj))
  x$Yr <- paste0(x$Yj, " / ", x$Yl)
  x$Xr <- paste0(x$Xj, " / ", x$Xl)
  oy <- order(match(x$Yl, y_parts), match(x$Yj, y_parts))
  ox <- order(match(x$Xl, x_parts), match(x$Xj, x_parts))
  x$Yr <- factor(x$Yr, levels = unique(x$Yr[oy]))
  x$Xr <- factor(x$Xr, levels = rev(unique(x$Xr[ox])))

  # lines between blocks with the same denominator
  block_lines <- function(ratio, den) {
    den <- den[match(levels(ratio), ratio)]
    which(den[-1] != den[-length(den)]) + 0.5
  }

  # blue diamonds < 0, grey dot at 0, red dots > 0
  bin_cols <- hcl.colors(n_bins, "Blue-Red 3")
  bin_cols[half + 1] <- "grey"
  bin_sizes <- sqrt(abs(seq(-half, half))) + 1
  bin_sizes <- bin_sizes / max(bin_sizes) * max_size
  bin_shapes <- rep(c(18, 20), c(half, half + 1))
  names(bin_cols) <- names(bin_sizes) <- names(bin_shapes) <- bin_labels
  # intToUtf8(177) is the plus-minus sign, kept ASCII for R CMD check
  legend_title <- paste0(
    "share ratio\nelasticity\n",
    intToUtf8(177),
    " ",
    signif(width / 2, 3)
  )

  ggplot2::ggplot(x) + # show.legend = TRUE: a legend key also for the bins without data
    ggplot2::geom_point(
      ggplot2::aes(x = Yr, y = Xr, col = SREB, pch = SREB, size = SREB),
      show.legend = TRUE
    ) +
    ggplot2::coord_fixed() +
    ggplot2::geom_hline(yintercept = block_lines(x$Xr, x$Xl)) +
    ggplot2::geom_vline(xintercept = block_lines(x$Yr, x$Yl)) +
    ggplot2::scale_color_manual(
      values = bin_cols,
      breaks = rev(bin_labels),
      drop = FALSE
    ) +
    ggplot2::scale_size_manual(
      values = bin_sizes,
      breaks = rev(bin_labels),
      drop = FALSE
    ) +
    ggplot2::scale_shape_manual(
      values = bin_shapes,
      breaks = rev(bin_labels),
      drop = FALSE
    ) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::labs(
      x = "share ratios of Y",
      y = "share ratios of X",
      col = legend_title,
      shape = legend_title,
      size = legend_title
    ) +
    ggplot2::theme_bw() +
    ggplot2::guides(color = ggplot2::guide_legend()) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = -75, hjust = 1, vjust = .5),
      legend.justification = c(0, 0)
    )
}

# columns of the plotted data used inside ggplot2::aes()
#' @importFrom utils globalVariables
NULL
globalVariables(c("Yr", "Xr", "SREB"))
