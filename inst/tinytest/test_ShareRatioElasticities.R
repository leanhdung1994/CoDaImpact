# ---- ShareRatioElasticities(): named numeric Xdir ----------------------------
sr_res <- lmCoDa(
  ilr(cbind(left, right, extreme_right)) ~
    ilr(cbind(Age_1839, Age_4064, Age_65plus)) + asset_owner_rate,
  data = election
)
sr_parm <- "cbind(Age_1839, Age_4064, Age_65plus)"
sr_a <- ShareRatioElasticities(
  sr_res,
  sr_parm,
  Xdir = c(Age_1839 = 1, Age_4064 = 2, Age_65plus = 3)
)
sr_b <- ShareRatioElasticities(
  sr_res,
  sr_parm,
  Xdir = c(Age_65plus = 3, Age_1839 = 1, Age_4064 = 2)
)
expect_equal(sr_a$SRE, sr_b$SRE)
expect_error(
  ShareRatioElasticities(sr_res, sr_parm, Xdir = c(a = 1, b = 2, c = 3)),
  pattern = "names of Xdir"
)

# ---- ShareRatioElasticities(): values of (35) --------------------------------
sr_B <- t(coef(sr_res, space = "clr", split = TRUE)[[2]]) # rows: parts of Y, columns: parts of X

# Xdir = NULL: (B[j, m] - B[l, m] - B[j, m'] + B[l, m']) / 2
sr_null <- ShareRatioElasticities(sr_res, sr_parm)
sr_exp <- mapply(
  function(yj, yl, xj, xl) {
    (sr_B[yj, xj] - sr_B[yl, xj] - sr_B[yj, xl] + sr_B[yl, xl]) / 2
  },
  sr_null$Yj,
  sr_null$Yl,
  sr_null$Xj,
  sr_null$Xl
)
expect_equal(unname(sr_null$SRE), unname(sr_exp))

# fixed direction u: sum_k (B[j, k] - B[l, k]) log u_k / log(u_m / u_m')
sr_u <- c(.2, .5, .3)
sr_fix <- ShareRatioElasticities(sr_res, sr_parm, Xdir = sr_u)
sr_lu <- setNames(log(sr_u), colnames(sr_B))
sr_exp2 <- mapply(
  function(yj, yl, xj, xl) {
    sum((sr_B[yj, ] - sr_B[yl, ]) * sr_lu) / (sr_lu[[xj]] - sr_lu[[xl]])
  },
  sr_fix$Yj,
  sr_fix$Yl,
  sr_fix$Xj,
  sr_fix$Xl
)
expect_equal(unname(sr_fix$SRE), unname(sr_exp2))
