#' A hazard generator for a discrete-time weibull distribution
#'
#' @param x A vector of positive real numbers at which to evaluate the hazard
#' @param shape The shape parameter
#' @param scale The scale parameter
#'
#' @return A vector of hazards evaluated at each element of argument 'x'
#' @export
#'
#' @examples
weihaz <- function(
  x,
  shape,
  scale) {

  # Check argument types
  checkmate::assert_numeric(shape)
  checkmate::assert_numeric(scale)
  checkmate::assert_numeric(x)

  # Defend against incorrect argument dimensions
  if (length(shape) != 1) stop("'shape' must be of length 1")
  if (length(scale) != 1) stop("'scale' must be of length 1")
  if (length(x) < 1) stop("'x' must be a vector of at least length 1")

  # Defend against incorrect argument magnitudes
  if (shape <= 0) stop("shape parameter must be positive")
  if (scale <= 0) stop("scale parameter must be positive")
  if (range(x)[1] <= 0) stop("x parameter must all be positive")

  return((shape/scale)*(x/scale)^(shape-1))
}
