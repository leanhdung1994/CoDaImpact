# Tests of plot.ShareRatioElasticities(), skipped without ggplot2
exit_if_not(requireNamespace("ggplot2", quietly = TRUE))

pl_res <- lmCoDa(
  ilr(cbind(left, right, extreme_right)) ~
    ilr(cbind(Age_1839, Age_4064, Age_65plus)),
  data = election
)
pl_parm <- "cbind(Age_1839, Age_4064, Age_65plus)"
pl_sre <- ShareRatioElasticities(pl_res, pl_parm)

# ---- class of the result ----------------------------------------------------
expect_inherits(pl_sre, "ShareRatioElasticities")
expect_inherits(pl_sre, "data.frame")

pl_p <- plot(pl_sre, n_bins = 11)
expect_inherits(pl_p, "ggplot")
expect_equal(nrow(pl_p$data), nrow(pl_sre))

# every bin has a legend key, also when it has no data
expect_true(isTRUE(pl_p$layers[[1]]$show.legend))

# ---- bins ---------------------------------------------------------------------
# default width: all bins exist, the largest |SRE| is in an outer bin
pl_bins <- pl_p$data$SREB
expect_equal(nlevels(pl_bins), 11)
pl_outer <- levels(pl_bins)[c(1, 11)]
expect_true(as.character(pl_bins[which.max(abs(pl_sre$SRE))]) %in% pl_outer)

# n_bins changes the number of bins
expect_equal(nlevels(plot(pl_sre, n_bins = 7)$data$SREB), 7)

# a tiny width puts every nonzero value in an outer bin
pl_tiny <- plot(pl_sre, n_bins = 11, width = 1e-8)$data$SREB
expect_true(all(as.character(pl_tiny) %in% levels(pl_tiny)[c(1, 11)]))

# the legend title shows width / 2
pl_w <- plot(pl_sre, n_bins = 11, width = 0.045)
pl_labs <- if (exists("get_labs", envir = asNamespace("ggplot2"))) {
  ggplot2::get_labs(pl_w) # ggplot2 >= 3.5.2
} else {
  pl_w$labels
}
expect_true(grepl("0.0225", pl_labs$colour, fixed = TRUE))

# max_size scales every symbol: max_size = 10 doubles the sizes of max_size = 5
expect_equal(
  ggplot2::layer_data(plot(pl_sre, n_bins = 11, max_size = 10))$size,
  2 * ggplot2::layer_data(plot(pl_sre, n_bins = 11))$size
)

# a value on the edge between two bins goes to the lower bin
pl_edge <- pl_sre
pl_edge$SRE <- rep_len(c(0.5, -0.5), nrow(pl_edge))
pl_edge_bins <- plot(pl_edge, width = 1, n_bins = 3)$data$SREB
expect_equal(
  as.character(pl_edge_bins),
  rep_len(c("0", "-1"), nrow(pl_edge))
)

# ---- non-finite values ----------------------------------------------------------
# towards a summit, the ratio of the two other parts does not change
pl_summit <- ShareRatioElasticities(pl_res, pl_parm, Xdir = "Age_65plus")
expect_message(plot(pl_summit, n_bins = 11), pattern = "not finite")
pl_ps <- suppressMessages(plot(pl_summit, n_bins = 11))
expect_equal(nrow(pl_ps$data), sum(is.finite(pl_summit$SRE)))

# ---- invalid arguments ----------------------------------------------------------
expect_error(plot(pl_sre), pattern = "n_bins must be given")
expect_error(plot(pl_sre, n_bins = 10), pattern = "n_bins")
expect_error(plot(pl_sre, n_bins = 1), pattern = "n_bins")
expect_error(plot(pl_sre, n_bins = 11, width = -1), pattern = "width")
expect_error(plot(pl_sre, n_bins = 11, max_size = -1), pattern = "max_size")
expect_error(
  plot(pl_sre[, c("Yj", "Yl", "Xj", "Xl")], n_bins = 11),
  pattern = "columns"
)
