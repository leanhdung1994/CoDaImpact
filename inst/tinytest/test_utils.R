# Tests of the internal helpers in R/utils.R

# ---- check_Xdir(): arguments -------------------------------------------------
ut_opts <- c("A", "B", "C")
expect_error(CoDaImpact:::check_Xdir("A", ut_opts, normalize = NA))
expect_error(CoDaImpact:::check_Xdir("A", ut_opts, normalize = c(TRUE, FALSE)))
expect_error(CoDaImpact:::check_Xdir("A", "A"))   # one part is not a composition
