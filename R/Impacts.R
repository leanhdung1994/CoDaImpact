#' @title Computation of elasticities in CoDa regression models
#'
#' @description
#' This function computes elasticities and semi-elasticities for assessing the
#' impact of a given covariate X on the dependent Y in a regression model. We
#' have to distinguish four cases:
#'  - Y and X are both compositional: this leads to an elasticity
#'  - Y is compositional and X is scalar: this leads to a semi-elasticity
#'  - Y is scalar and X is compositional: this leads to a semi-elasticity
#'  - Y and X are both scalar: this case is not implemented, as it is the same as in ordinary linear regression, yielding a constant marginal effect
#'
#' @details
#' The mathematical foundation for elasticity computations in CoDa models comes
#' from Morais and Thomas-Agnan (2021). Dargel and Thomas-Agnan (2024) present
#' further results and illustrations.
#'
#' @param object an object of class "lmCoDa"
#' @param Xvar a character indicating the name of one explanatory variable
#' @param obs a numeric that refers to the indicator of one observation
#'
#' @return a matrix
#' @author
#'   - Lukas Dargel
#'   - Rodrigue Nasr
#' @export
#' @references
#'   - Dargel, Lukas and Christine Thomas-Agnan, “Pairwise share ratio interpretations of compositional regression models”, Computational Statistics & Data Analysis 195 (2024), p. 107945
#'   - Morais, Joanna and Christine Thomas-Agnan. "Impact of covariates in compositional models and simplicial derivatives." Austrian Journal of Statistics 50.2 (2021): 1-15.
#' @examples
#' res <- lmCoDa(
#'   YIELD ~ PRECIPITATION + ilr(TEMPERATURES),
#'   data = head(rice_yields, 20)
#' )
#' Elasticities(res, Xvar = "TEMPERATURES")
#'
Elasticities <- function(object, Xvar = NULL, obs = 1) {
  stopifnot(
    is(object, "lmCoDa"),
    is.character(Xvar) && length(Xvar) == 1,
    is.numeric(obs) && isTRUE(obs >= 1) && obs <= nobs(object)
  )

  trSry <- object$trSry
  Xvar <- check_Xvar(Xvar, trSry, return_type = "pos")
  Xcoef <- trSry$COEF_CLR[[Xvar]]

  check <- "Elasticities are only meaningful if X or Y are compositional!"
  YX_is_compo <- c("" != trSry$LR_TRAN[c(1, Xvar)], use.names = FALSE)
  if (identical(YX_is_compo, c(FALSE, FALSE))) {
    stop(check)
  }

  if (identical(YX_is_compo, c(FALSE, TRUE))) {
    return(Xcoef)
  }

  # when Y is compositional
  Dy <- trSry$D[[1]]
  Wz <- diag(Dy) -
    as(fitted(object, space = "simplex")[rep(obs, Dy), ], "matrix")
  imp <- Xcoef %*% t(Wz)
  colnames(imp) <- colnames(Xcoef)
  attr(imp, "obs") <- obs
  return(imp)
}

#' @title Impacts in CoDa regression models through pairwise share ratios
#'
#' @description
#' This function computes the left-hand sides of equations (29), (32) and (35)
#' of Dargel and Thomas-Agnan (2024), with the first equality of each:
#' - Y is scalar and X is compositional, (29): the change of Y per relative
#' change of the share ratio Xj / Xl, along the direction u
#' - Y is compositional and X is scalar, (32): the relative change of the share
#' ratio Yj / Yl per unit change of X
#' - Y and X are both compositional, (35): the relative change of the share
#' ratio Yj / Yl per relative change of the share ratio Xj / Xl, along the
#' direction u
#'
#' @param object an object of class "lmCoDa"
#' @param Xvar a character indicating the name of one explanatory variable
#' @param Xdir a numeric vector, a single character, or `NULL`, only used when X
#'   is compositional:
#'   - if numeric, `Xdir` is the direction u in the simplex
#'   - if character, `Xdir` is one part of X and u is the direction towards it
#'   - if `NULL`, u is the direction that multiplies Xj by e and divides Xl by e,
#'     for each pair (Xj, Xl)
#'
#' @return
#' a data.frame with one row per pair of distinct parts and the impact in the
#' column `FI`
#' @export
#' @references
#' - Dargel, Lukas and Christine Thomas-Agnan, "Pairwise share ratio interpretations of compositional regression models", Computational Statistics & Data Analysis 195 (2024), p. 107945
Impacts <- function(object, Xvar, Xdir = NULL) {
  stopifnot(
    is(object, "lmCoDa"),
    is.character(Xvar) && length(Xvar) == 1,
    is.null(Xdir) || is.character(Xdir) || is.numeric(Xdir)
  )

  trSry <- object$trSry
  Xpos <- check_Xvar(Xvar, trSry, return_type = "pos")
  y_is_compo <- trSry$LR_TRAN[[1]] != ""
  x_is_compo <- trSry$LR_TRAN[[Xpos]] != ""

  if (!y_is_compo && !x_is_compo) {
    stop("Impacts are only meaningful if X or Y are compositional!")
  }

  if (!x_is_compo && !is.null(Xdir)) {
    stop("Xdir must be NULL when X is scalar!")
  }

  # clr coefficients: one row per part of X, one column per part of Y
  coef_clr <- trSry$COEF_CLR[[Xpos]]
  x_names <- rownames(coef_clr)
  y_names <- colnames(coef_clr)

  # direction u for the pair (Xj, Xl) of parts of X
  get_u <- function(Xj, Xl) {
    if (!is.null(Xdir)) {
      return(check_Xdir(Xdir, Xopts = x_names, normalize = FALSE))
    }
    u <- rep(1, length(x_names))
    u[Xj] <- exp(1)
    u[Xl] <- exp(-1)
    u / sum(u)
  }

  # pairs of distinct parts of X and of Y
  pairs_x <- expand.grid(Xj = seq_along(x_names), Xl = seq_along(x_names))
  pairs_x <- pairs_x[pairs_x$Xj != pairs_x$Xl, ]
  pairs_y <- expand.grid(Yj = seq_along(y_names), Yl = seq_along(y_names))
  pairs_y <- pairs_y[pairs_y$Yj != pairs_y$Yl, ]

  # (29): scalar Y, compositional X
  if (!y_is_compo) {
    clr_beta <- coef_clr[, 1]
    result <- pairs_x
    result$FI <- NA_real_
    for (r in seq_len(nrow(result))) {
      Xj <- result$Xj[r]
      Xl <- result$Xl[r]
      u <- get_u(Xj, Xl)
      clr_u <- unclass(clr(u))
      result$FI[r] <- sum(clr_beta * clr_u) / log(u[Xj] / u[Xl])
    }
    result$Xj <- x_names[result$Xj]
    result$Xl <- x_names[result$Xl]
    rownames(result) <- NULL
    return(result)
  }

  # (32): compositional Y, scalar X
  if (!x_is_compo) {
    gamma <- unclass(trSry$COEF_SIMPLEX[[Xpos]])[1, ]
    result <- pairs_y
    result$FI <- log(gamma[result$Yj] / gamma[result$Yl])
    result$Yj <- y_names[result$Yj]
    result$Yl <- y_names[result$Yl]
    rownames(result) <- NULL
    return(result)
  }

  # (35): compositional Y and X
  result <- merge(pairs_y, pairs_x)
  result$FI <- NA_real_
  for (r in seq_len(nrow(result))) {
    Yj <- result$Yj[r]
    Yl <- result$Yl[r]
    Xj <- result$Xj[r]
    Xl <- result$Xl[r]
    u <- get_u(Xj, Xl)
    numerator <- sum((coef_clr[, Yj] - coef_clr[, Yl]) * log(u))
    result$FI[r] <- numerator / log(u[Xj] / u[Xl])
  }
  result$Yj <- y_names[result$Yj]
  result$Yl <- y_names[result$Yl]
  result$Xj <- x_names[result$Xj]
  result$Xl <- x_names[result$Xl]
  return(result)
}
