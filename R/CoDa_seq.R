#' A sequence connecting two points in the simplex by a simplex line
#'
#' @details
#' The sequence is evenly spaced and corresponds to a straight line in the
#' simplex geometry. If no end point is provided the line will connect the
#' initial point with the first summit of the simplex. Since exact zeros are not
#' handled by the ilr they are replaced by a small constant.
#'
#' @param comp_from A numeric vector, representing the initial composition
#' @param comp_to A numeric vector, representing the final composition
#' @param n_steps An integer, indicating the number of steps used to go from
#'   comp_from to comp_to
#' @param add_opposite A logical, if `TRUE` the path in the opposite direction
#'   is added
#'
#' @importFrom compositions acomp norm
#' @importFrom methods is as
#' @return
#' A data.frame where each row corresponds to one compositional vector, closed
#' to 1
#' @seealso simplex_increment
#' @export
#'
#' @author Lukas Dargel
#' @examples
#' # path to the first summit of the simplex
#' start_comp <- c(A = .4, B = .35, C = .25)
#' compositions::plot.acomp(CoDa_seq(start_comp))
#' compositions::plot.acomp(CoDa_seq(start_comp, add_opposite = TRUE))
#'
#' # path to an edge of the simplex
#' end_comp <- c(0, .8, .2)
#' compositions::plot.acomp(CoDa_seq(start_comp, end_comp))
#' compositions::plot.acomp(CoDa_seq(start_comp, end_comp, add_opposite = TRUE))
CoDa_seq <- function(comp_from, comp_to, n_steps = 100, add_opposite = FALSE) {
  epsilon <- 1e-10
  if (missing(comp_to)) {
    comp_to <- rep(epsilon, length(comp_from))
    names(comp_to) <- names(comp_from)
    comp_to[1] <- 1
  }

  stopifnot(
    is.numeric(comp_from),
    is.numeric(comp_to),
    length(comp_from) >= 2,
    length(comp_from) == length(comp_to),
    all(comp_from >= 0),
    all(comp_to >= 0),
    sum(comp_from) > 0,
    sum(comp_to) > 0,
    is.null(names(comp_to)) || identical(names(comp_to), names(comp_from)),
    is.numeric(n_steps),
    length(n_steps) == 1,
    round(n_steps) >= 1,
    isTRUE(add_opposite) || isFALSE(add_opposite)
  )

  n_steps <- round(n_steps)

  check <- " must be a single compositional vector!"
  if (isTRUE(nrow(comp_from) > 1)) {
    stop("comp_from", check)
  }
  if (isTRUE(nrow(comp_to) > 1)) {
    stop("comp_to", check)
  }
  comp_from <- as(comp_from, "vector")
  comp_to <- as(comp_to, "vector")

  # close both compositions, and replace zeros
  comp_from <- comp_from / sum(comp_from)
  comp_to <- comp_to / sum(comp_to)

  zero_values_present <- any(comp_from == 0, comp_to == 0)
  if (zero_values_present) {
    warn <- "Zero values in `comp_from` or `comp_to` are not permitted and were replaced by %s (after closing to 1)!"
    warning(sprintf(warn, epsilon))
    comp_from[comp_from == 0] <- epsilon
    comp_from <- comp_from / sum(comp_from)
    comp_to[comp_to == 0] <- epsilon
    comp_to <- comp_to / sum(comp_to)
  }

  x <- acomp(comp_from)
  y <- acomp(comp_to)
  comp_dir <- y - x
  k <- seq(-n_steps * add_opposite, n_steps)
  result <- lapply(k / n_steps, function(s) unclass(x + s * comp_dir))
  result <- do.call("rbind", result)
  row.names(result) <- k

  if (!is.null(names(comp_from))) {
    colnames(result) <- names(comp_from)
  }

  Anorm_dir <- norm(comp_dir)
  if (Anorm_dir < epsilon) {
    stop("`comp_from` and `comp_to` must be different compositions!")
  }

  result <- as.data.frame(result)
  attr(result, "step_size") <- Anorm_dir / n_steps
  attr(result, "direction") <- (1 / Anorm_dir) * comp_dir
  return(result)
}

#' Create a linear path in the simplex by defining a direction and a step size
#'
#' @details
#' The function is very similar to [CoDa_seq()]. However, instead of drawing a
#' line between a starting point and an end point, it uses only a starting point
#' and a direction.
#'
#' @param comp_direc A numeric vector, defining a direction in the simplex
#' @param comp_from A numeric vector, an initial point in the simplex - defaults
#'   to the origin of the simplex (the closure of the vector of ones)
#' @param step_size A numeric, indicating the step size
#' @param n_steps A numeric, indicating the number of steps to be taken from
#'   `comp_from`
#' @param add_opposite A logical, if `TRUE` steps in the opposite direction are
#'   also computed
#' @param dir_from_start A logical, if `TRUE` the direction is calculated from
#'   the difference between `comp_from` and `comp_direc`
#'
#' @return
#' A data.frame where each row corresponds to one compositional vector, closed
#' to 1
#' @importFrom compositions acomp norm
#' @seealso CoDa_seq
#' @export
#'
#' @author Lukas Dargel
#' @examples
#' # three steps that go from the origin towards the defined direction
#' comp_direc <- c(A = .4, B = .35, C = .25)
#' CoDa_path(comp_direc, n_steps = 3)
#'
#' # we can draw the path that is defined by this direction
#' comp_direc <- c(A = .4, B = .35, C = .25)
#' compositions::plot.acomp(CoDa_path(comp_direc, n_steps = 10))
#' compositions::plot.acomp(CoDa_path(comp_direc, n_steps = 100))
#' compositions::plot.acomp(CoDa_path(comp_direc, add_opposite = TRUE))
#'
#' # using the same direction we can draw a new path that does not go through the origin
#' comp_direc <- c(A = .4, B = .35, C = .25)
#' comp_from <- c(.7, .2, .1)
#' compositions::plot.acomp(CoDa_path(comp_direc, comp_from, n_steps = 10))
#' compositions::plot.acomp(CoDa_path(comp_direc, comp_from, n_steps = 100))
#' compositions::plot.acomp(CoDa_path(
#'   comp_direc,
#'   comp_from,
#'   add_opposite = TRUE
#' ))
#'
#' # the balanced composition does not define a direction by itself
#' comp_origin <- c(A = 1 / 3, B = 1 / 3, C = 1 / 3) # corresponds to a zero vector in real space
#' try(CoDa_path(comp_origin, comp_from, add_opposite = TRUE))
#'
#' # with the dir_from_start option the direction is derived
#' # from the simplex line connecting two compositions
#' path_origin <- CoDa_path(
#'   comp_direc = comp_origin,
#'   comp_from = comp_from,
#'   add_opposite = TRUE,
#'   dir_from_start = TRUE,
#'   step_size = .1
#' )
#' compositions::plot.acomp(path_origin)
#' compositions::plot.acomp(comp_origin, add = TRUE, col = "blue", pch = 19)
#' compositions::plot.acomp(comp_from, add = TRUE, col = "red", pch = 19)
#'
CoDa_path <- function(
  comp_direc,
  comp_from,
  step_size = 0.01,
  n_steps = 100,
  add_opposite = FALSE,
  dir_from_start = FALSE
) {
  D <- length(comp_direc)
  if (missing(comp_from)) {
    comp_from <- rep(1 / D, D)
    names(comp_from) <- names(comp_direc)
  }

  stopifnot(
    is.numeric(comp_direc),
    is.numeric(comp_from),
    length(comp_direc) >= 2,
    length(comp_from) == length(comp_direc),
    all(comp_direc >= 0),
    all(comp_from >= 0),
    sum(comp_direc) > 0,
    sum(comp_from) > 0,
    is.numeric(step_size),
    length(step_size) == 1,
    is.finite(step_size),
    step_size > 0,
    is.numeric(n_steps),
    length(n_steps) == 1,
    round(n_steps) >= 1,
    isTRUE(add_opposite) || isFALSE(add_opposite),
    isTRUE(dir_from_start) || isFALSE(dir_from_start)
  )

  n_steps <- round(n_steps)

  if (isTRUE(nrow(comp_from) > 1)) {
    stop("comp_from must be a single compositional vector!")
  }
  if (isTRUE(nrow(comp_direc) > 1)) {
    stop("comp_direc must be a single compositional vector!")
  }
  if (is.matrix(comp_from)) {
    comp_from <- comp_from[1, ]
  }
  if (is.matrix(comp_direc)) {
    comp_direc <- comp_direc[1, ]
  }

  if (is.null(names(comp_from))) {
    names(comp_from) <- names(comp_direc)
  }
  stopifnot(
    is.null(names(comp_direc)) || identical(names(comp_from), names(comp_direc))
  )

  # close both compositions, and replace zeros
  epsilon <- 1e-10
  comp_from <- comp_from / sum(comp_from)
  comp_direc <- comp_direc / sum(comp_direc)

  zero_values_present <- any(comp_from == 0, comp_direc == 0)
  if (zero_values_present) {
    warn <- "Zero values in `comp_from` or `comp_direc` are not permitted and were replaced by %s (after closing to 1)!"
    warning(sprintf(warn, epsilon))
    comp_from[comp_from == 0] <- epsilon
    comp_from <- comp_from / sum(comp_from)
    comp_direc[comp_direc == 0] <- epsilon
    comp_direc <- comp_direc / sum(comp_direc)
  }

  x <- acomp(comp_from)
  comp_dir <- acomp(comp_direc)
  if (dir_from_start) {
    comp_dir <- comp_dir - x
  }

  Anorm_dir <- norm(comp_dir)
  if (Anorm_dir < epsilon) {
    stop("The simplex direction is invalid!")
  }
  comp_dir <- (1 / Anorm_dir) * comp_dir

  # calculate the end point and use CoDa_seq
  comp_to <- x + (n_steps * step_size) * comp_dir
  result <- CoDa_seq(comp_from, unclass(comp_to), n_steps, add_opposite)
  return(result)
}
