#' Compute share ratio elasticities for CoDa models
#'
#' In CoDa models with compositional dependent variable (Y), share ratio
#' elasticities (SRE) allow to interpret the influence of compositional
#' explanatory variables (X). The interpretation is analogous to usual
#' elasticities:
#'   - When the share ratio of X increases by 1% in a given direction of the simplex specified by `Xdir`, each share ratio of Y increases by its SRE%
#'   - The main difference to usual elasticities is that, since X is compositional
#'     the change of X must be specified in terms of a direction in the simplex.
#'
#' For a scalar Y or a scalar X, see [Impacts()].
#'
#'
#' @details
#' The function computes the left-hand side of equation (35) of Dargel and
#' Thomas-Agnan (2024) for every combination of a pair of distinct parts
#' \eqn{(Y_j, Y_l)}{(Yj, Yl)} of \eqn{Y} and a pair of distinct parts
#' \eqn{(X_j, X_l)}{(Xj, Xl)} of \eqn{X}, the two pairs being chosen
#' independently. The result has one row per combination, given in the columns
#' `Yj`, `Yl`, `Xj` and `Xl`, and the share ratio elasticity in the column
#' `SRE`. The direction \eqn{u} is:
#' - the same for all rows when `Xdir` is given
#' - dependent on the pair \eqn{(X_j, X_l)}{(Xj, Xl)} when `Xdir` is `NULL`
#'
#' More details on this interpretation can be found in Dargel and Thomas-Agnan
#' (2024) and in the accompanying vignette.
#'
#' @inheritParams VariationScenario
#' @param Xdir a numeric vector, the name of one part of X, or `NULL`:
#'   - if numeric, `Xdir` is a fixed direction in the simplex
#'   - if the name of a part, the direction is towards the summit of this part
#'   - if `NULL`, the direction depends on the pair (Xj, Xl): it multiplies Xj
#'     by e and divides Xl by e
#' @return a data.frame
#'
#' @author Lukas Dargel
#' @references
#'   - Dargel, Lukas and Christine Thomas-Agnan, “Pairwise share ratio interpretations of compositional regression models”, Computational Statistics & Data Analysis 195 (2024), p. 107945
#' @importFrom compositions acomp
#' @export
#' @examples
#' ### regression model with a compositional dependent and at least one compositional covariate
#' res <- lmCoDa(
#'   ilr(cbind(left, right, extreme_right)) ~
#'     ilr(cbind(Educ_BeforeHighschool, Educ_Highschool, Educ_Higher)),
#'   data = head(election, 20)
#' )
#' educ_comp <- "cbind(Educ_BeforeHighschool, Educ_Highschool, Educ_Higher)"
#'
#' ## case 1: direction towards the summit "Educ_Higher"
#' SRE1 <- ShareRatioElasticities(res, Xvar = educ_comp, Xdir = "Educ_Higher")
#'
#' # when Educ_Higher / Educ_BeforeHighschool increases by 1% in this direction,
#' # right / left increases by about SRE%
#' subset(SRE1, Yj == "right" & Yl == "left" &
#'   Xj == "Educ_Higher" & Xl == "Educ_BeforeHighschool")
#'
#' # in this direction Educ_Highschool / Educ_BeforeHighschool does not change,
#' # so the SRE is not defined
#' subset(SRE1, Yj == "right" & Yl == "left" &
#'   Xj == "Educ_Highschool" & Xl == "Educ_BeforeHighschool")
#'
#' ## case 2: the same direction as a numeric vector
#' SRE2 <- ShareRatioElasticities(res, Xvar = educ_comp, Xdir = exp(c(0, 0, 1)))
#' identical(SRE1, SRE2)
#'
#' ## case 3: Xdir = NULL, for each pair (Xj, Xl) the direction multiplies Xj by e
#' ## and divides Xl by e
#' SRE3 <- ShareRatioElasticities(res, Xvar = educ_comp, Xdir = NULL)
#' subset(SRE3, Yj == "right" & Yl == "left" &
#'   Xj == "Educ_Highschool" & Xl == "Educ_BeforeHighschool")
ShareRatioElasticities <- function(object, Xvar, Xdir = NULL) {
  stopifnot(
    is(object, "lmCoDa"),
    is.character(Xvar) && length(Xvar) == 1,
    is.null(Xdir) || is.character(Xdir) || is.numeric(Xdir)
  )

  trSry <- object$trSry
  Xpos <- check_Xvar(Xvar, trSry, "pos")

  check <- "Share ratio elasticities are only meaningful when X and Y are compositional!"
  Dx <- trSry$D[[Xpos]]
  Dy <- trSry$D[[1]]
  if (Dx == 0 || Dy == 0) {
    stop(check)
  }

  # allocate for results in vectorized form
  resX <- expand.grid(Xj = seq(Dx), Xl = seq(Dx))
  resY <- expand.grid(Yj = seq(Dy), Yl = seq(Dy))
  resXY <- merge(resY[resY$Yj != resY$Yl, ], resX[resX$Xj != resX$Xl, ])
  resXY <- cbind(resXY, "SRE" = NA_real_)

  y_names <- colnames(trSry$COEF_CLR[[Xpos]])
  x_names <- rownames(trSry$COEF_CLR[[Xpos]])

  if (!is.null(Xdir)) {
    Xdir <- check_Xdir(Xdir = Xdir, Xopts = x_names, normalize = FALSE)
  }

  # direction u in the simplex for the pair (Xj, Xl) of parts of X
  get_u <- function(Xj, Xl) {
    if (!is.null(Xdir)) {
      return(Xdir)
    }
    # u^(Xj) minus u^(Xl): multiplies Xj by e and divides Xl by e
    u <- acomp(exp(seq_len(Dx) == Xj)) - acomp(exp(seq_len(Dx) == Xl))
    return(unclass(u))
  }

  B <- t(trSry$COEF_CLR[[Xpos]])
  for (r in seq_len(nrow(resXY))) {
    Yj <- resXY$Yj[r]
    Yl <- resXY$Yl[r]
    Xj <- resXY$Xj[r]
    Xl <- resXY$Xl[r]
    u <- get_u(Xj, Xl)

    # share ratio elasticity (35) for the direction u
    numerator <- sum((B[Yj, ] - B[Yl, ]) * log(u))
    denominator <- log(u[Xj] / u[Xl])
    resXY$SRE[r] <- numerator / denominator
  }

  resXY$Yj <- y_names[resXY$Yj]
  resXY$Yl <- y_names[resXY$Yl]
  resXY$Xj <- x_names[resXY$Xj]
  resXY$Xl <- x_names[resXY$Xl]
  return(resXY)
}
