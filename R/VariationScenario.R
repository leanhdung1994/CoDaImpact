#' @title Scenarios for variation in CoDa regressions models
#'
#' @description
#' Scenarios of this type are illustrated in Dargel and Thomas-Agnan (2024).
#' They allow to evaluate how the response variable (Y) in a CoDa model would
#' evolve under a hypothetical scenario for linear changes in one explanatory
#' variable (X). When the changing explanatory variable is compositional the
#' term "linear" is understood with respect to the geometry of the simplex.
#'
#' @details
#' The scenario of X consists of the points X(h) for the steps h = k * inc_size,
#' with k from -n_steps (or 0) to n_steps. For a scalar X, X(h) = X0 + h. For a
#' compositional X, X(h) is the straight line in the simplex that starts at X0
#' and goes in the direction Xdir. The corresponding values of Y are computed
#' with the prediction formula, where we exploit the fact that only one variable
#' changes.
#'
#'
#' @param object an object of class "lmCoDa"
#' @param Xvar a character indicating the name of the explanatory variable that
#'   changes
#' @param Xdir either character or numeric, to indicate the direction in which
#'   Xvar should change
#'   - when character this should be one of the components of X, in which case
#'     the direction points towards the corresponding vertex of the simplex
#'   - when numeric this should be a positive vector, which is used as a
#'     direction in the simplex
#'   - (when Xvar refers to a scalar variable this argument is ignored)
#' @param obs a positive integer indicating the row of the observation used for
#'   the scenario
#' @param inc_size a positive numeric indicating the step size of the scenario
#'   of X
#' @param n_steps a positive integer indicating the number of steps in the
#'   scenario, in each direction
#' @param add_opposite a logical, if `TRUE` the scenario also includes changes
#'   in the opposite direction
#' @param normalize_Xdir a logical, if `TRUE` the direction `Xdir` is scaled to
#'   have an Aitchison norm of 1, allowing to interpret `inc_size` as the
#'   Aitchison distance
#' @return
#' a data.frame containing the scenario of X and the corresponding predicted
#' values of Y
#' @author Lukas Dargel
#' @references
#'   - Dargel, Lukas and Christine Thomas-Agnan, “Pairwise share ratio interpretations of compositional regression models”, Computational Statistics & Data Analysis 195 (2024), p. 107945
#' @importFrom compositions acomp
#' @export
#' @examples
#' # ---- model with scalar response ----
#' res <- lmCoDa(
#'   YIELD ~ PRECIPITATION + ilr(TEMPERATURES),
#'   data = head(rice_yields, 20)
#' )
#' VariationScenario(res, Xvar = "TEMPERATURES", Xdir = "MEDIUM", n_steps = 5)
#' VariationScenario(res, Xvar = "PRECIPITATION", n_steps = 5)
#'
#' # ---- model with compositional response ----
#' res <- lmCoDa(
#'   ilr(cbind(left, right, extreme_right)) ~
#'     ilr(cbind(Age_1839, Age_4064)) +
#'       ilr(cbind(Educ_BeforeHighschool, Educ_Highschool, Educ_Higher)) +
#'       log(unemp_rate),
#'   data = head(election)
#' )
#'
#' VariationScenario(
#'   res,
#'   Xvar = "cbind(Age_1839,Age_4064)",
#'   Xdir = "Age_1839",
#'   n_steps = 5
#' )
#' VariationScenario(res, "log(unemp_rate)", n_steps = 5)
#'
VariationScenario <- function(
  object,
  Xvar,
  Xdir,
  obs = 1,
  inc_size = .1,
  n_steps = 100,
  add_opposite = TRUE,
  normalize_Xdir = TRUE
) {
  stopifnot(
    is(object, "lmCoDa"),
    is.character(Xvar) && length(Xvar) == 1,
    missing(Xdir) || is.character(Xdir) || is.numeric(Xdir),
    is.numeric(obs) && isTRUE(obs >= 1) && obs <= nobs(object),
    obs == round(obs),
    is.numeric(inc_size) && length(inc_size) == 1,
    is.finite(inc_size) && inc_size > 0,
    is.numeric(n_steps) && length(n_steps) == 1,
    is.finite(n_steps) && n_steps >= 1,
    n_steps == round(n_steps),
    isTRUE(add_opposite) || isFALSE(add_opposite),
    isTRUE(normalize_Xdir) || isFALSE(normalize_Xdir)
  )

  clo2 <- function(x) {
    x / rowSums(x)
  }
  trSry <- object$trSry

  # get X0
  Xvar <- check_Xvar(Xvar, trSry, "NAME_COORD")
  Dx <- trSry$D[[Xvar]]
  scalar_x <- Dx == 0
  Kx <- trSry$LR_BASE_K[[Xvar]]
  X0 <- object$model[[Xvar]]
  X0 <- if (scalar_x) X0[obs] else t(X0[obs, ]) %*% t(Kx) # work with clr in compositional case

  # get Y0
  scalar_y <- trSry$D[[1]] == 0
  Ky <- trSry$LR_BASE_K[[1]]
  Y0 <- fitted(object)
  Y0 <- if (scalar_y) Y0[obs] else Y0[obs, ] %*% t(Ky) # work with clr in compositional case

  inc_seq <- seq(from = -n_steps * add_opposite, to = n_steps)
  if (scalar_x) {
    # IDEA allow for "log-linear" sequences
    Xcoef <- trSry$COEF_CLR[[Xvar]]
    Xscenario <- inc_seq * inc_size
    Yscenario <- if (scalar_y) {
      Y0
    } else {
      Y0[rep(1, length(Xscenario)), , drop = FALSE]
    }
    Yscenario <- Xscenario %*% Xcoef + Yscenario
    Xscenario <- X0 + Xscenario
  }

  if (!scalar_x) {
    if (missing(Xdir)) {
      stop("Xdir must be supplied when Xvar is compositional!")
    }
    Xdir <- check_Xdir(Xdir, colnames(X0), normalize_Xdir)
    Xcoef <- trSry$COEF_CLR[[Xvar]]

    # scenario of X in the simplex: x0 + h * u, with h = k * inc_size
    x0 <- clrInv(X0[1, ])
    u <- acomp(Xdir)
    h <- inc_seq * inc_size
    Xscenario <- lapply(h, function(s) unclass(x0 + s * u))
    Xscenario <- do.call("rbind", Xscenario)
    colnames(Xscenario) <- colnames(X0)

    # response in clr (or scalar response): Y0 + h * clr(u) B
    Yscenario <- if (scalar_y) Y0 else Y0[rep(1, length(h)), ]
    Yscenario <- h %*% (unclass(clr(u)) %*% Xcoef) + Yscenario
  }

  if (!scalar_x) {
    Xvar <- "X"
  }
  if (!scalar_y) {
    Yvar <- "Y"
  } else {
    Yvar <- rownames(trSry)[1]
  }
  if (!scalar_y) {
    Yscenario <- clo2(exp(Yscenario))
  }
  result <- data.frame(row.names = as.character(inc_seq))
  result[[Yvar]] <- Yscenario
  result[[Xvar]] <- Xscenario
  attr(result, "inc_size") <- inc_size
  attr(result, "obs") <- obs
  return(result)
}
