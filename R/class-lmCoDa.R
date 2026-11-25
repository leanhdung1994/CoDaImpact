# ---- constructors -----------------------------------------------------------

#' Estimating CoDa regression models
#'
#' This is a thin wrapper around [lm()] followed by [ToSimplex()], which allows
#' to create a lmCoDa object in one step.
#'
#' @param formula as in [lm()]
#' @param data as in [lm()]
#' @param ... arguments passed on to [lm()]
#' @return
#' an object of class "lm" and "lmCoDa" if the formula includes at least one
#' log-transformation
#'
#' @author Lukas Dargel
#' @seealso [lm()], [ToSimplex()], [compositions::ilr()], [compositions::alr()]
#' @export
#' @examples
#' # regression model with a compositional dependent and at least one compositional covariate
#' res <- lmCoDa(
#'   ilr(cbind(left, right, extreme_right)) ~
#'     ilr(cbind(Educ_BeforeHighschool, Educ_Highschool, Educ_Higher)),
#'   data = head(election, 20)
#' )
#'
#' # regression model with a scalar dependent and at least one compositional covariate
#' res <- lmCoDa(
#'   YIELD ~ PRECIPITATION + ilr(TEMPERATURES),
#'   data = head(rice_yields, 20)
#' )
#'
lmCoDa <- function(formula, data, ...) {
  cl <- match.call()
  cl[[1]] <- as.symbol("lm")
  res <- eval(cl)
  cl[[1]] <- as.symbol("lmCoDa")
  res$call <- cl
  return(ToSimplex(res))
}

#' Converting Linear Models to CoDa models
#'
#' @description
#' The function converts the output of a "lm" to the "lmCoDa" class, which
#' offers additional tools for the interpretation of a CoDa regression model.
#' Most of the work is done by the [transformationSummary()] function, which has
#' its own documentation page, but should be reserved for internal use.
#'
#'
#' @inherit lmCoDa return examples
#' @inheritParams VariationScenario
#' @seealso [lm()], [lmCoDa()]
#' @author
#'   - Lukas Dargel
#'   - Rodrigue Nasr
#' @export
#' @examples
#' # regression model with a compositional dependent and at least one compositional covariate
#' res <- lm(
#'   ilr(cbind(left, right, extreme_right)) ~
#'     ilr(cbind(Educ_BeforeHighschool, Educ_Highschool, Educ_Higher)),
#'   data = head(election, 20)
#' )
#' res <- ToSimplex(res)
#'
#' # regression model with a scalar dependent and at least one compositional covariate
#' res <- lm(
#'   YIELD ~ PRECIPITATION + ilr(TEMPERATURES),
#'   data = head(rice_yields, 20)
#' )
#' res <- ToSimplex(res)
ToSimplex <- function(object) {
  trSry <- transformationSummary(object)
  if (all(trSry[, 'LR_TRAN'] == "")) {
    return(object)
  }

  object$trSry <- trSry
  class(object) <- c('lmCoDa', class(object))
  if (trSry[1, 'LR_TRAN'] == "") {
    return(object)
  }

  invTran <- match.fun(paste0(trSry[1, 'LR_TRAN'], "Inv"))
  attributes(object$fitted.values)$orig <- invTran(object$fitted.values)
  attributes(object$residuals)$orig <- invTran(object$residuals)
  return(object)
}

# ---- methods ----------------------------------------------------------------
#' Coefficients of CoDa regression models
#'
#' Works as [coef()] for lm objects. The `space` argument returns the
#' coefficients in clr space or in the simplex instead of the log-ratio
#' coordinates used in the model.
#'
#' @param object an object of class "lmCoDa"
#' @param space `NULL` (default) for the coefficients in the log-ratio
#'   coordinates used in the model, or one of `"clr"` and `"simplex"`
#' @param split logical, if `TRUE` the coefficients are reported as a list
#'   instead of a matrix, where list structure reflects the explanatory
#'   variables of the model
#' @param ... not used
#' @return
#' the coefficients of the explanatory variables. Each scalar variable has one
#' row; a compositional variable with l parts has l - 1 rows with `space = NULL`
#' and l rows with `space = "clr"` or `"simplex"`.
#'
#' @author Lukas Dargel
#' @exportS3Method
coef.lmCoDa <- function(object, space = NULL, split = FALSE, ...) {
  stopifnot(is.null(space) || space %in% c("clr", "simplex"))
  type <- if (is.null(space)) {
    "COEF_COORD"
  } else {
    c("clr" = "COEF_CLR", "simplex" = "COEF_SIMPLEX")[space]
  }
  cfs <- object$trSry[[type]][-1]
  cfs <- if (split) cfs else Reduce("rbind", cfs)
  return(cfs)
}

#' Fitted values of CoDa regression models
#'
#' Works as [fitted()] for lm objects. When the response Y is compositional, the
#' `space` argument returns the fitted values in clr space or in the simplex
#' instead of the log-ratio coordinates used in the model.
#'
#' @inheritParams predict.lmCoDa
#' @param ... not used
#' @return the fitted values, in the same form as for [predict.lmCoDa()]
#' @importFrom stats fitted
#' @author Lukas Dargel
#' @exportS3Method
fitted.lmCoDa <- function(object, space = NULL, ...) {
  stopifnot(is.null(space) || space %in% c("clr", "simplex"))
  fity <- object$fitted.values
  if (is.null(space)) {
    return(fity)
  }

  Ky <- object$trSry[["LR_BASE_K"]][[1]]
  if (length(Ky) == 0) {
    stop("The space argument can only be used for Y compositional models!")
  }

  if (space == "simplex") {
    return(attr(fity, "orig"))
  }

  if (space == "clr") return(fity %*% t(Ky))
}

#' Predictions from CoDa regression models
#'
#' Works as [predict.lm()]. When the response Y is compositional, the `space`
#' argument returns the predictions in clr space or in the simplex instead of
#' the log-ratio coordinates used in the model.
#'
#' @param object an object of class "lmCoDa"
#' @param space `NULL` (default) for the log-ratio coordinates used in the
#'   model, or one of `"clr"` and `"simplex"`. Only available when Y is
#'   compositional.
#' @param ... passed on to [predict.lm()], for example `newdata`
#' @return
#' - with `space = NULL`, the same as [predict.lm()]: a vector for a scalar Y,
#'   a matrix of log-ratio coordinates for a compositional Y
#' - with `space = "clr"`, a matrix of clr coordinates with one row per
#'   observation and one column per part of Y
#' - with `space = "simplex"`, a matrix of compositions with one row per
#'   observation, closed to 1
#' @author Lukas Dargel
#' @exportS3Method
predict.lmCoDa <- function(object, space = NULL, ...) {
  stopifnot(is.null(space) || space %in% c("clr", "simplex"))
  pred <- NextMethod("predict")
  if (is.null(space)) {
    return(pred)
  }

  Ky <- object$trSry[["LR_BASE_K"]][[1]]
  if (length(Ky) == 0) {
    stop("The space argument can only be used for Y compositional models!")
  }

  if (space == "simplex") {
    return(clrInv(pred %*% t(Ky)))
  }

  if (space == "clr") return(pred %*% t(Ky))
}

#' Residuals of CoDa regression models
#'
#' Works as [residuals()] for lm objects. When the response Y is compositional,
#' the `space` argument returns the residuals in clr space or in the simplex
#' instead of the log-ratio coordinates used in the model.
#'
#' @inheritParams predict.lmCoDa
#' @param ... not used
#' @return
#' the residuals, in the same form as for [predict.lmCoDa()]. In the simplex,
#' the residual of an observation is the composition Y minus the fitted
#' composition (the perturbation difference).
#' @author Lukas Dargel
#' @exportS3Method
residuals.lmCoDa <- function(object, space = NULL, ...) {
  stopifnot(is.null(space) || space %in% c("clr", "simplex"))
  resi <- object$residuals
  if (is.null(space)) {
    return(resi)
  }

  Ky <- object$trSry[["LR_BASE_K"]][[1]]
  if (length(Ky) == 0) {
    stop("The space argument can only be used for Y compositional models!")
  }

  if (space == "simplex") {
    return(attr(resi, "orig"))
  }

  if (space == "clr") return(resi %*% t(Ky))
}
