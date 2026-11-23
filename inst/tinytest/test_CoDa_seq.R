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

# ---- CoDa_seq() and CoDa_path(): paths in the simplex ------------------------

# evenly spaced straight line of proportions
cs_x <- c(A = 2, B = 3, C = 5)
cs_y <- c(A = .5, B = .3, C = .2)
cs_res <- CoDa_seq(cs_x, cs_y, n_steps = 4, add_opposite = TRUE)

expect_equal(rownames(cs_res), as.character(-4:4))
expect_equal(unname(rowSums(cs_res)), rep(1, 9))
expect_equal(unlist(cs_res["0", ]), cs_x / sum(cs_x))
expect_equal(unlist(cs_res["4", ]), cs_y / sum(cs_y))

# equal steps in clr coordinates = straight, evenly spaced line
cs_clr <- unclass(compositions::clr(as.matrix(cs_res)))
cs_steps <- diff(cs_clr)
expect_equal(
  cs_steps,
  matrix(cs_steps[1, ], nrow(cs_steps), ncol(cs_steps), byrow = TRUE),
  check.attributes = FALSE
)
expect_equal(attr(cs_res, "step_size"), sqrt(sum(cs_steps[1, ]^2)))
expect_equal(as.numeric(compositions::norm(attr(cs_res, "direction"))), 1)

# two equal compositions, zeros
expect_error(CoDa_seq(c(.2, .3, .5), c(2, 3, 5)), pattern = "different compositions")
expect_warning(
  CoDa_seq(c(A = .4, B = .35, C = .25), c(A = 0, B = .8, C = .2)),
  pattern = "after closing to 1"
)

# CoDa_path(): steps of the given Aitchison length
cs_path <- CoDa_path(c(A = .4, B = .35, C = .25), step_size = .1, n_steps = 5)
expect_equal(attr(cs_path, "step_size"), .1)
expect_equal(unname(rowSums(cs_path)), rep(1, 6))
expect_error(CoDa_path(c(A = 1 / 3, B = 1 / 3, C = 1 / 3)), pattern = "invalid")
