#' @title Effects of infinitesimal changes in CoDa models
#'
#' @description
#' This function allows to evaluate how a change in an explanatory variable
#' impacts the response variable in a CoDa regression model. The changes are
#' calculated with the first-order approximation presented in Dargel and
#' Thomas-Agnan (2024). Changes in the response variable are provided as a
#' data.frame and the underlying changes in the explanatory variable are given
#' as attributes.
#'
#'
#' @inheritParams VariationScenario
#' @param inc_rate a numeric between 0 and 1, only used when `Xdir` is the name
#'   of a part of X. It is the relative change of this part, and replaces
#'   `inc_size` by the step size that gives this change.
#' @param Ytotal a positive numeric indicating the total of Y, used to convert
#'   the variation of the parts into units (only used when Y is compositional)
#' @return
#' a data.frame with one column per part of Y (one column for a scalar Y) and
#' the rows
#' - the initial value of Y at the observation `obs`
#' - the new value of Y after the step, to first order
#' - the elasticity (or semi-elasticity) of Y with respect to X along `Xdir`
#' - the variation of Y in %, in % points (compositional Y) and in units
#' (scalar Y, or compositional Y with `Ytotal` different from 1)
#'
#' The attributes `X(0)` and `X(h)` give the value of X before and after the
#' step, `Xdir` the direction, `inc_size` the step size, and `inc_rates` the
#' relative change of each part of X. `Xdir` and `inc_rates` are `NULL` for a
#' scalar X.
#'
#' @author
#'   - Lukas Dargel
#'   - Rodrigue Nasr
#' @references
#'   - Dargel, Lukas and Christine Thomas-Agnan, “Pairwise share ratio interpretations of compositional regression models”, Computational Statistics & Data Analysis 195 (2024), p. 107945
#' @export
#' @examples
#' # XY-compositional model
#' res <- lmCoDa(
#'   ilr(cbind(left, right, extreme_right)) ~
#'     ilr(cbind(Educ_BeforeHighschool, Educ_Highschool, Educ_Higher)),
#'   data = head(election, 20)
#' )
#'
#' # Focus on changes in the education composition
#' educ_comp <- "cbind(Educ_BeforeHighschool, Educ_Highschool, Educ_Higher)"
#'
#' # ... changes towards a summit (higher share of people with lower education)
#' VariationTable(res, educ_comp, Xdir = "Educ_BeforeHighschool")
#'
#' # ... same changes using a compositional vector as direction
#' VariationTable(res, educ_comp, Xdir = c(.5, .25, .25))
#'
#' # ... changes in a more general direction and for a different observation
#' VariationTable(res, educ_comp, Xdir = c(.35, .45, .10), obs = 2)
#'
VariationTable <- function(
  object,
  Xvar,
  Xdir,
  obs = 1,
  inc_size = .1,
  inc_rate = NULL,
  Ytotal = 1,
  normalize_Xdir = TRUE
) {
  stopifnot(
    is(object, "lmCoDa"),
    is.character(Xvar) && length(Xvar) == 1,
    missing(Xdir) || is.character(Xdir) || is.numeric(Xdir),
    is.numeric(obs) && isTRUE(obs >= 1) && obs <= nobs(object),
    is.numeric(inc_size) && length(inc_size) == 1,
    is.null(inc_rate) || (0 < inc_rate && inc_rate < 1),
    is.numeric(Ytotal) && length(Ytotal) == 1
  )

  trSry <- object$trSry
  if (all(trSry$LR_TRAN[c(1, Xvar)] == "")) {
    stop("Variation tables are only meaningful if X or Y are compositional!")
  }

  # get X0
  Xvar <- check_Xvar(Xvar, trSry, "NAME_COORD")
  Dx <- trSry$D[[Xvar]]
  scalar_x <- Dx == 0
  X0 <- object$model[[Xvar]]
  X0 <- if (scalar_x) X0[obs] else attr(X0, "orig")[obs, ]

  # get Y0
  scalar_y <- trSry$D[[1]] == 0
  Y0 <- if (scalar_y) {
    structure(fitted(object)[obs], names = names(object$model)[1])
  } else {
    as(fitted(object, space = "simplex")[obs, ], "vector")
  }

  # define elasticities
  elasti <- Elasticities(object, trSry$NAME_SIMPLEX[[Xvar]], obs)

  if (scalar_x & !is.null(inc_rate)) {
    warning("The argument inc_rate is ignored when X is a classical variable!")
  }
  if (!scalar_x) {
    # for compositional X we need to account for the direction in which X changes
    vertex_dir <- is.character(Xdir)
    Xdir <- check_Xdir(Xdir, names(X0), normalize_Xdir)
    elasti <- log(Xdir) %*% elasti

    # link between alpha (ink_rate) and h (ink_size)
    inc_rates <- log(Xdir) - sum(as(X0, "vector") * log(Xdir))
    if (vertex_dir && !is.null(inc_rate)) {
      inc_size <- inc_rate / inc_rates[which.max(Xdir)]
    }
    inc_rates <- inc_rates * inc_size
  }

  Ytotal <- if (scalar_y) 1 else Ytotal
  Y1 <- Y0 * (1 + elasti * inc_size)
  YD <- Y1 - Y0
  explanation <- c(
    "Y0" = if (scalar_y) "Initial value" else "Initial parts",
    "Y1" = if (scalar_y) "New value" else "New parts",
    "elasti" = if (scalar_y || scalar_x) "Semi elasticity" else "Elasticity",
    "YD/Y0*100" = "Variation in %",
    "YD*100" = if (scalar_y) NULL else "Variation in % points",
    "YD*Ytotal" = if (scalar_y || Ytotal != 1) "Variation in units" else NULL
  )

  evalString <- function(s) {
    eval(str2lang(s), environment())
  }
  explanation <- explanation[!is.null(explanation)]
  result <- lapply(names(explanation), evalString)
  result <- data.frame(
    Reduce("rbind", result),
    row.names = explanation,
    check.names = FALSE
  )

  if (scalar_x) {
    Xdir <- inc_rates <- NULL
  }
  attr(result, 'X(0)') <- X0
  attr(result, 'X(h)') <- if (scalar_x) {
    X0 + inc_size
  } else {
    as(X0, "vector") * (1 + inc_rates)
  }
  attr(result, 'Xdir') <- Xdir
  attr(result, 'inc_size') <- inc_size
  attr(result, 'inc_rates') <- inc_rates
  return(result)
}
