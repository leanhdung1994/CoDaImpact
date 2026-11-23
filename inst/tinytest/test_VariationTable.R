# ---- Y - compositional case ------------------------------------------------
res <- lmCoDa(
  ilr(cbind(left, right, extreme_right)) ~
    ilr(cbind(Age_1839, Age_4064)) +
    ilr(cbind(Educ_BeforeHighschool, Educ_Highschool, Educ_Higher)) +
    log(unemp_rate),
  data = head(election, 20))

# compositional X
expect_equal(
  VariationTable(res, Xvar = "cbind(Age_1839, Age_4064)", Xdir = "Age_1839"),
  VariationTable(res, Xvar = "cbind(Age_1839, Age_4064)", Xdir = c(2,1)),
  info = "Coherence between preferential and general direction.",
  check.attributes = FALSE)

resX <- VariationTable(res, Xvar = "cbind(Educ_BeforeHighschool, Educ_Highschool, Educ_Higher)", Xdir = c(.5,.4,.1))
expect_true("elasticity" %in% tolower(rownames(resX)))
expect_true(sum(attr(resX, "X(0)")) == 1)
expect_true(sum(attr(resX, "X(h)")) == 1)
expect_true(sum(attr(resX, "Xdir")) == 1)
expect_true(sum(resX["Initial parts",]) == 1)
expect_equivalent(sum(resX["New parts",]),1)
expect_equivalent(sum(resX["Variation in % points",]),0)
rm(resX)

# scalar X
resZ <- VariationTable(res, Xvar = "log(unemp_rate)")
expect_true("semi elasticity" %in% tolower(rownames(resZ)))
expect_equivalent(sum(resZ["Initial parts",]), 1)
expect_equivalent(sum(resZ["New parts",]), 1)
expect_equivalent(sum(resZ["Variation in % points",]),0)
rm(resZ, res)


# ---- Y scalar ---------------------------------------------------------------
res <- lmCoDa(YIELD ~ PRECIPITATION + ilr(TEMPERATURES), data = head(rice_yields, 20))

# compositional X
resX <- VariationTable(res, Xvar = "TEMPERATURES", Xdir = "HIGH")
expect_true(colnames(resX) == "YIELD")
expect_true("semi elasticity" %in% tolower(rownames(resX)))
expect_true(sum(attr(resX, "X(0)")) == 1)
expect_true(sum(attr(resX, "Xdir")) == 1)
expect_true(resX["Initial value",] == res$fitted.values[[1]])
expect_true(isTRUE(abs(resX["Variation in units",]) > 0))
rm(resX)


# scalar X
expect_error(VariationTable(res, Xvar = "PRECIPITATION"),pattern =  "meaningful")

# ---- argument checks ---------------------------------------------------------
vt_fit <- lmCoDa(
  ilr(cbind(left, right, extreme_right)) ~
    ilr(cbind(Age_1839, Age_4064)) + log(unemp_rate),
  data = head(election, 20)
)
vt_X <- "cbind(Age_1839, Age_4064)"
expect_error(VariationTable(vt_fit, vt_X, "Age_1839", obs = 1.5))
expect_error(VariationTable(vt_fit, vt_X, "Age_1839", inc_size = 0))
expect_error(VariationTable(vt_fit, vt_X, "Age_1839", inc_size = NA_real_))
expect_error(VariationTable(vt_fit, vt_X, "Age_1839", inc_size = Inf))
expect_error(VariationTable(vt_fit, vt_X, "Age_1839", inc_rate = c(.1, .2)))
expect_error(VariationTable(vt_fit, vt_X, "Age_1839", Ytotal = NA_real_))
expect_error(VariationTable(vt_fit, vt_X, "Age_1839", Ytotal = -1))
expect_error(VariationTable(vt_fit, vt_X, "Age_1839", normalize_Xdir = NA))
expect_error(VariationTable(vt_fit, vt_X), pattern = "Xdir must be supplied")
expect_warning(
  VariationTable(vt_fit, vt_X, Xdir = c(2, 1), inc_rate = .05),
  pattern = "inc_rate is ignored when Xdir is numeric"
)

# ---- scalar X and scalar Y: no variation table --------------------------------
vt_fit_s <- lmCoDa(
  YIELD ~ PRECIPITATION + ilr(TEMPERATURES),
  data = head(rice_yields, 20)
)
expect_error(
  VariationTable(vt_fit_s, Xvar = "PRECIPITATION"),
  pattern = "Variation tables are only meaningful"
)

# ---- first-order values, compared with VariationScenario() --------------------
vt_u <- c(.2, .5, .3)

# scalar Y: the model is linear in h, so the new value is exact
vt_h <- .5
vt_tab <- VariationTable(vt_fit_s, "TEMPERATURES", Xdir = vt_u, obs = 2,
  inc_size = vt_h)
vt_vs <- VariationScenario(vt_fit_s, "TEMPERATURES", Xdir = vt_u, obs = 2,
  inc_size = vt_h, n_steps = 1, add_opposite = FALSE)
vt_y <- as.numeric(vt_vs[[1]])
expect_equal(vt_tab["Initial value", 1], vt_y[1])
expect_equal(vt_tab["New value", 1], vt_y[2], info = "additive for scalar Y")
expect_equal(vt_tab["Variation in units", 1], vt_y[2] - vt_y[1])
expect_equal(
  vt_tab["New value", 1],
  vt_tab["Initial value", 1] + vt_h * vt_tab["Semi elasticity", 1]
)

# compositional Y: agrees with the exact scenario to first order
vt_h <- 1e-4
vt_tab <- VariationTable(vt_fit, vt_X, Xdir = "Age_1839", inc_size = vt_h,
  Ytotal = 1000)
vt_vs <- VariationScenario(vt_fit, vt_X, Xdir = "Age_1839", inc_size = vt_h,
  n_steps = 1, add_opposite = FALSE)
vt_y <- unclass(vt_vs[[1]])
expect_equal(as.numeric(vt_tab["Initial parts", ]), as.numeric(vt_y[1, ]))
expect_equal(as.numeric(vt_tab["New parts", ]), as.numeric(vt_y[2, ]),
  tolerance = 1e-6)
expect_equal(
  as.numeric(vt_tab["Variation in units", ]),
  1000 * as.numeric(vt_tab["New parts", ] - vt_tab["Initial parts", ])
)
expect_equal(rownames(vt_tab), c("Initial parts", "New parts", "Elasticity",
  "Variation in %", "Variation in % points", "Variation in units"))
expect_equal(colnames(vt_tab), c("left", "right", "extreme_right"))
