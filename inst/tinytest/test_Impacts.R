# ---- Y - compositional ------------------------------------------------------
res <- lmCoDa(
  ilr(cbind(left, right, extreme_right)) ~
    ilr(cbind(Age_1839, Age_4064)) +
      ilr(cbind(Educ_BeforeHighschool, Educ_Highschool, Educ_Higher)) +
      unemp_rate,
  data = election[1:20, ]
)

expect_equal(
  Elasticities(res, Xvar = "cbind(Age_1839, Age_4064)", obs = 10) * 0,
  {
    res0 <- matrix(
      0,
      2,
      3,
      dimnames = list(
        c("Age_1839", "Age_4064"),
        c("left", "right", "extreme_right")
      )
    )
    attr(res0, "obs") <- 10
    res0
  }
)

expect_equal(Elasticities(res, Xvar = "unemp_rate", obs = 10) * 0, {
  res0 <- matrix(
    0,
    1,
    3,
    dimnames = list(c("unemp_rate"), c("left", "right", "extreme_right"))
  )
  attr(res0, "obs") <- 10
  res0
})

# ---- Y - scalar -------------------------------------------------------------
res <- lmCoDa(
  YIELD ~ PRECIPITATION + ilr(TEMPERATURES),
  data = rice_yields[1:20, ]
)

expect_error(
  Elasticities(res, Xvar = "PRECIPITATION"),
  pattern = "Elasticities"
)

expect_equal(Elasticities(res, Xvar = "TEMPERATURES", obs = 10) * 0, {
  res0 <- matrix(
    0,
    3,
    1,
    dimnames = list(c("LOW", "MEDIUM", "HIGH"), c("YIELD"))
  )
  res0
})

# a vector of names is rejected
expect_error(Elasticities(res, Xvar = c("TEMPERATURES", "PRECIPITATION")))

# ---- Impacts(): pairwise share ratio impacts ---------------------------------

# (29) scalar Y, compositional X: b' clr(u) / log(u_j / u_l), with b = Elasticities()
im_res <- lmCoDa(
  YIELD ~ PRECIPITATION + ilr(TEMPERATURES),
  data = rice_yields[1:20, ]
)
im_u <- c(LOW = .2, MEDIUM = .5, HIGH = .3)
im_SE <- Elasticities(im_res, Xvar = "TEMPERATURES")[, 1]
im_FI <- Impacts(im_res, Xvar = "TEMPERATURES", Xdir = im_u)
expect_equal(nrow(im_FI), 6)
for (r in seq_len(nrow(im_FI))) {
  expect_equal(
    unname(im_FI$FI[r]),
    sum(im_SE * log(im_u)) / log(im_u[[im_FI$Xj[r]]] / im_u[[im_FI$Xl[r]]])
  )
}

# (32) compositional Y, scalar X: difference of two semi-elasticities
ic_res <- lmCoDa(
  ilr(cbind(left, right, extreme_right)) ~
    ilr(cbind(Age_1839, Age_4064, Age_65plus)) + asset_owner_rate,
  data = election
)
ic_FI <- Impacts(ic_res, Xvar = "asset_owner_rate")
ic_SE <- Elasticities(ic_res, Xvar = "asset_owner_rate", obs = 1)[1, ]
expect_equal(nrow(ic_FI), 6)
expect_equal(unname(ic_FI$FI), unname(ic_SE[ic_FI$Yj] - ic_SE[ic_FI$Yl]))

# (35) compositional Y and X, default direction: equals ShareRatioElasticities()
ic_parm <- "cbind(Age_1839, Age_4064, Age_65plus)"
ic_FI35 <- Impacts(ic_res, Xvar = ic_parm)
ic_SRE <- ShareRatioElasticities(ic_res, Xvar = ic_parm)
expect_equal(nrow(ic_FI35), 36)
expect_equal(unname(ic_FI35$FI), unname(ic_SRE$SRE))

# errors
expect_error(
  Impacts(im_res, Xvar = "PRECIPITATION"),
  pattern = "only meaningful"
)
expect_error(
  Impacts(ic_res, Xvar = "asset_owner_rate", Xdir = "left"),
  pattern = "Xdir"
)
