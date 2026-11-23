
# ---- ShareRatioElasticities(): named numeric Xdir ----------------------------
sr_res <- lmCoDa(
  ilr(cbind(left, right, extreme_right)) ~
    ilr(cbind(Age_1839, Age_4064, Age_65plus)) + asset_owner_rate,
  data = election
)
sr_parm <- "cbind(Age_1839, Age_4064, Age_65plus)"
sr_a <- ShareRatioElasticities(sr_res, sr_parm, Xdir = c(Age_1839 = 1, Age_4064 = 2, Age_65plus = 3))
sr_b <- ShareRatioElasticities(sr_res, sr_parm, Xdir = c(Age_65plus = 3, Age_1839 = 1, Age_4064 = 2))
expect_equal(sr_a$SRE, sr_b$SRE)
expect_error(
  ShareRatioElasticities(sr_res, sr_parm, Xdir = c(a = 1, b = 2, c = 3)),
  pattern = "names of Xdir"
)
