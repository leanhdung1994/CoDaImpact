# ---- test CoDa_seq ------------------------------------------------------
expect_true({
  res <- c(A =.4,B = .3, C= .3)
  res <- CoDa_seq(res,n_steps = 2)
  is.data.frame(res)
})

expect_true({
  res <- c(A =.4,B = .3, C= .3)
  res <- CoDa_seq(res,n_steps = 2)
  all(colnames(res) == c("A","B","C"))
})

expect_true({
  res <- c(A =.4,B = .3, C= .3)
  res <- CoDa_seq(res,n_steps = 2,add_opposite = TRUE)
  all(dim(res) == c(5,3))
})


expect_true({
  res <- c(A =.4,B = .3, C= .3)
  res <- CoDa_seq(res,n_steps = 2)
  all(abs(res["0",] - t(c(A =.4,B = .3, C= .3))) < .Machine$double.eps)
})


expect_true({
  res <- c(A =.4,B = .3, C= .3)
  res <- CoDa_seq(res,n_steps = 1,add_opposite = TRUE)
  all(row.names(res) == c("-1","0","1"))
})


# ---- test CoDa_path --------------------------------------------------
expect_true({
  res <- c(A =.4,B = .3, C= .3)
  res <- CoDa_path(res,n_steps = 1)
  is.data.frame(res)
})

expect_true({
  res <- c(A =.4,B = .3, C= .3)
  res <- CoDa_path(res,n_steps = 1)
  all(colnames(res) == c("A","B","C"))
})

expect_true({
  res <- c(A =.4,B = .3, C= .3)
  res <- CoDa_path(res,n_steps = 1,add_opposite = TRUE)
  all(dim(res) == c(3,3))
})

expect_true({
  res <- c(A =.4,B = .3, C= .3)
  res <- CoDa_path(res,n_steps = 1,add_opposite = TRUE)
  all(res[2,] == 1/c(3,3))
})


expect_true({
  res <- c(A =.4,B = .3, C= .3)
  res <- CoDa_path(res,n_steps = 1,add_opposite = FALSE,step_size = 1)
  abs(1 - attr(res, "step_size")) < 1e12
})


expect_equal({
  vec <- c(A = 2, B = 1, C = 1)
  CoDa_path(vec,n_steps = 1,add_opposite = TRUE)
},
{
  vec <- c(A = 9, B = 1, C = 1)
  CoDa_path(vec,n_steps = 1,add_opposite = TRUE)
}, info = "normalization of directions to unit length")


# ---- CoDa_seq() and CoDa_path(): argument checks -----------------------------

cs_x <- c(A = .2, B = .3, C = .5)
cs_y <- c(A = .5, B = .3, C = .2)
expect_error(CoDa_seq(cs_x, cs_y, add_opposite = NA))
expect_error(CoDa_seq(cs_x, cs_y, add_opposite = c(TRUE, FALSE)))
expect_error(CoDa_seq(cs_x, cs_y, n_steps = c(1, 2)))
expect_error(CoDa_seq(cs_x, cs_y, n_steps = 0))
expect_error(CoDa_seq(c(A = -.2, B = .7, C = .5), cs_y))
expect_error(CoDa_seq(c(A = NA, B = .5, C = .5), cs_y))
expect_error(CoDa_seq(c(0, 0, 0), c(.5, .3, .2)))
expect_error(CoDa_seq(c("a", "b"), c("c", "d")))

cs_z <- c(A = .4, B = .35, C = .25)
expect_error(CoDa_path(c("a", "b")))
expect_error(CoDa_path(c(0, 0, 0)))
expect_error(CoDa_path(cs_z, step_size = 0))
expect_error(CoDa_path(cs_z, step_size = -1))
expect_error(CoDa_path(cs_z, step_size = Inf))
expect_error(CoDa_path(cs_z, add_opposite = NA))
expect_error(CoDa_path(cs_z, dir_from_start = NA))
expect_error(CoDa_path(cs_z, c(B = .2, A = .3, C = .5)))

# default end point of CoDa_seq() is (1, eps, eps), whatever comp_from is
cs_res <- CoDa_seq(c(A = .2, B = .3, C = .5), n_steps = 2)
cs_end <- cs_res[nrow(cs_res), ]
expect_equal(cs_end$B / cs_end$C, 1, tolerance = 1e-6)
