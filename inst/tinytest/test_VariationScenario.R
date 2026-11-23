expect_equal({
  res <- lmCoDa(
    ilr(cbind(left, right, extreme_right)) ~
      ilr(cbind(Age_1839, Age_4064)) +
      ilr(cbind(Educ_BeforeHighschool, Educ_Highschool, Educ_Higher)) +
      log(unemp_rate),
    data = head(election,20))

  list(
    VariationScenario(res, "cbind(Age_1839, Age_4064)",Xdir = "Age_1839", n_steps = 5, add_opposite = FALSE),
    VariationScenario(res, "log(unemp_rate)", n_steps = 5, add_opposite = FALSE))
},
{
  res <- lmCoDa(
    alr(cbind(left, right, extreme_right)) ~
      alr(cbind(Age_1839, Age_4064)) +
      alr(cbind(Educ_BeforeHighschool, Educ_Highschool, Educ_Higher)) +
      log(unemp_rate),
    data = head(election,20))
  list(
    VariationScenario(res, "cbind(Age_1839, Age_4064)",Xdir = "Age_1839", add_opposite = FALSE, n_steps = 5),
    VariationScenario(res, "log(unemp_rate)",Xdir = "Age_1839",add_opposite = FALSE,n_steps = 5))
},info = "Y compo & (X compo + X scalar)")

expect_equal({
  res <- lmCoDa(YIELD ~ PRECIPITATION + ilr(TEMPERATURES), data = head(rice_yields,20))
  list(
    VariationScenario(res, "TEMPERATURES", Xdir = "MEDIUM", n_steps = 5, add_opposite = FALSE),
    VariationScenario(res, "PRECIPITATION", n_steps = 5, add_opposite = FALSE))
},
{
  res <- lmCoDa(YIELD ~ PRECIPITATION + alr(TEMPERATURES), data = head(rice_yields,20))
  list(
    VariationScenario(res, "TEMPERATURES", Xdir = "MEDIUM", n_steps = 5, add_opposite = FALSE),
    VariationScenario(res, "PRECIPITATION", n_steps = 5, add_opposite = FALSE))
},info = "Y scalar & (X compo + X scalar)")


# ---- argument checks ---------------------------------------------------------
vs_fit <- lmCoDa(
  YIELD ~ PRECIPITATION + ilr(TEMPERATURES),
  data = head(rice_yields, 20)
)
expect_error(VariationScenario(vs_fit, "TEMPERATURES", "MEDIUM", obs = 1.5))
expect_error(VariationScenario(vs_fit, "TEMPERATURES", "MEDIUM", inc_size = 0))
expect_error(VariationScenario(vs_fit, "TEMPERATURES", "MEDIUM", inc_size = -1))
expect_error(VariationScenario(vs_fit, "TEMPERATURES", "MEDIUM", inc_size = Inf))
expect_error(VariationScenario(vs_fit, "TEMPERATURES", "MEDIUM", inc_size = NA_real_))
expect_error(VariationScenario(vs_fit, "TEMPERATURES", "MEDIUM", n_steps = 0))
expect_error(VariationScenario(vs_fit, "TEMPERATURES", "MEDIUM", n_steps = 2.5))
expect_error(VariationScenario(vs_fit, "PRECIPITATION", normalize_Xdir = NA))
expect_error(
  VariationScenario(vs_fit, "TEMPERATURES"),
  pattern = "Xdir must be supplied"
)

# ---- compositional scenario: x0 + h * u in the simplex -----------------------
# (lmCoDa() looks up `data` from the package, so we use `election` directly)
vs_age <- c("Age_1839", "Age_4064", "Age_65plus")
vs_X <- "cbind(Age_1839, Age_4064, Age_65plus)"
vs_fit2 <- lmCoDa(
  ilr(cbind(left, right, extreme_right)) ~
    ilr(cbind(Age_1839, Age_4064, Age_65plus)) + log(unemp_rate),
  data = election
)

vs_clr <- function(x) log(x) - mean(log(x))
vs_clrInv <- function(z) exp(z) / sum(exp(z))
vs_u <- c(.2, .5, .3)
vs_h <- (-2:2) * .5
vs_x0 <- unlist(election[3, vs_age])
vs_tr <- vs_fit2$trSry
vs_y0 <- drop(unname(fitted(vs_fit2)[3, ] %*% t(vs_tr$LR_BASE_K[[1]])))  # clr
vs_B <- vs_tr$COEF_CLR[[which(vs_tr$NAME_SIMPLEX == vs_X)]]

# expected scenario of X and Y for a direction with clr coordinates d
vs_expect <- function(d) {
  list(
    X = t(sapply(vs_h, function(h) vs_clrInv(vs_clr(vs_x0) + h * d))),
    Y = t(sapply(vs_h, function(h) vs_clrInv(vs_y0 + h * (d %*% vs_B))))
  )
}

# normalize_Xdir = FALSE uses u as given (it was always normalised before)
vs_out <- VariationScenario(vs_fit2, vs_X, Xdir = vs_u, obs = 3,
  inc_size = .5, n_steps = 2, normalize_Xdir = FALSE)
vs_exp <- vs_expect(vs_clr(vs_u))
expect_equal(unname(unclass(vs_out[[2]])), unname(vs_exp$X), info = "X, raw Xdir")
expect_equal(unname(unclass(vs_out[[1]])), unname(vs_exp$Y), info = "Y, raw Xdir")

# normalize_Xdir = TRUE (default) uses u / ||u||_A: same values as before
vs_out <- VariationScenario(vs_fit2, vs_X, Xdir = vs_u, obs = 3,
  inc_size = .5, n_steps = 2)
vs_exp <- vs_expect(vs_clr(vs_u) / sqrt(sum(vs_clr(vs_u)^2)))
expect_equal(unname(unclass(vs_out[[2]])), unname(vs_exp$X), info = "X, unit Xdir")
expect_equal(unname(unclass(vs_out[[1]])), unname(vs_exp$Y), info = "Y, unit Xdir")
expect_equal(colnames(vs_out[[2]]), vs_age)

# ---- output: columns Y and X, attributes and print method --------------------
vs_out <- VariationScenario(vs_fit, "PRECIPITATION", n_steps = 2)
expect_equal(names(vs_out), c("Y", "X"))
expect_true(inherits(vs_out, "VariationScenario"))
expect_true(is.numeric(vs_out$Y) && is.null(dim(vs_out$Y)))
expect_equal(attr(vs_out, "Yvar"), "YIELD")
expect_equal(attr(vs_out, "Xvar"), "PRECIPITATION")
expect_stdout(print(vs_out), pattern = "scalar response `YIELD`")

vs_out <- VariationScenario(vs_fit2, vs_X, Xdir = "Age_1839", n_steps = 2)
expect_equal(names(vs_out), c("Y", "X"))
expect_true(is.matrix(vs_out$Y) && is.matrix(vs_out$X))
expect_equal(unname(rowSums(vs_out$Y)), rep(1, 5))
expect_equal(attr(vs_out, "Yvar"), "cbind(left, right, extreme_right)")
expect_equal(attr(vs_out, "Xvar"), vs_X)
expect_stdout(print(vs_out), pattern = "compositional explanatory variable `cbind(Age_1839", fixed = TRUE)

# a scalar regressor named Y no longer overwrites the predicted response
vs_fit3 <- lmCoDa(
  ilr(cbind(left, right, extreme_right)) ~ ilr(cbind(Age_1839, Age_4064)) + Y,
  data = transform(election, Y = unemp_rate)
)
vs_out <- VariationScenario(vs_fit3, "Y", n_steps = 2)
expect_equal(names(vs_out), c("Y", "X"))
expect_true(is.matrix(vs_out$Y))
expect_equal(vs_out$X, election$unemp_rate[1] + (-2:2) * .1)
