# Tests of the internal helpers in R/utils.R

# ---- check_Xdir(): arguments -------------------------------------------------
ut_opts <- c("A", "B", "C")
expect_error(CoDaImpact:::check_Xdir("A", ut_opts, normalize = NA))
expect_error(CoDaImpact:::check_Xdir("A", ut_opts, normalize = c(TRUE, FALSE)))
expect_error(CoDaImpact:::check_Xdir("A", "A")) # one part is not a composition

# ---- check_Xdir(): numeric direction -----------------------------------------
expect_error(CoDaImpact:::check_Xdir(c(0, .5, .5), ut_opts)) # zero part
expect_error(CoDaImpact:::check_Xdir(c(1, 1, 1), ut_opts))   # neutral element
expect_error(CoDaImpact:::check_Xdir(c(NA, .5, .5), ut_opts))
expect_error(CoDaImpact:::check_Xdir(c(.5, .5), ut_opts)) # wrong length

ut_u <- c(.2, .5, .3)
expect_equal(
  as.numeric(CoDaImpact:::check_Xdir(ut_u, ut_opts)),
  ut_u / sum(ut_u)
)

# a named numeric direction is matched to the parts by name, not by position
expect_equal(
  as.numeric(CoDaImpact:::check_Xdir(c(C = .3, A = .2, B = .5), ut_opts)),
  as.numeric(CoDaImpact:::check_Xdir(ut_u, ut_opts))
)
expect_error(
  CoDaImpact:::check_Xdir(c(a = .2, b = .5, c = .3), ut_opts),
  pattern = "names of Xdir"
)
expect_error(
  CoDaImpact:::check_Xdir(c(A = .2, .5, .3), ut_opts), # partly named
  pattern = "names of Xdir"
)

ut_n <- as.numeric(CoDaImpact:::check_Xdir(ut_u, ut_opts, normalize = TRUE))
expect_equal(sum(ut_n), 1)
expect_equal(as.numeric(compositions::norm(compositions::acomp(ut_n))), 1)
# same result as the previous code: ilr, divide by its norm, back
ut_ilr <- compositions::ilr(ut_u)
expect_equal(
  ut_n,
  as.numeric(unclass(compositions::ilrInv(ut_ilr / sqrt(sum(ut_ilr^2)))))
)

# ---- check_Xdir(): direction towards a summit --------------------------------
expect_error(CoDaImpact:::check_Xdir(c("A", "B"), ut_opts))
ut_s <- as.numeric(CoDaImpact:::check_Xdir("B", ut_opts))
expect_equal(ut_s, exp(c(0, 1, 0)) / sum(exp(c(0, 1, 0))))

ut_sn <- as.numeric(CoDaImpact:::check_Xdir("B", ut_opts, normalize = TRUE))
expect_equal(as.numeric(compositions::norm(compositions::acomp(ut_sn))), 1)
# same result as the previous code: e replaced by e^sqrt(l / (l - 1))
ut_e <- exp(c(0, 1, 0) * sqrt(3 / 2))
expect_equal(ut_sn, ut_e / sum(ut_e))

# ---- check_Xvar() ------------------------------------------------------------
ut_res <- lmCoDa(
  YIELD ~ PRECIPITATION + ilr(TEMPERATURES),
  data = rice_yields[1:20, ]
)
expect_error(CoDaImpact:::check_Xvar(NA_character_, ut_res$trSry))
expect_error(
  CoDaImpact:::check_Xvar("nonexistent", ut_res$trSry),
  pattern = "must be one of"
)
expect_equal(
  CoDaImpact:::check_Xvar("TEMPERATURES", ut_res$trSry, "NAME_COORD"),
  "ilr(TEMPERATURES)"
)
expect_equal(CoDaImpact:::check_Xvar(" TEMPERATURES ", ut_res$trSry, "pos"), 4L)
