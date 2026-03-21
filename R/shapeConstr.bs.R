################################################################################
#
# Shape constraint matrix method:
# bs method
#
################################################################################

#' @export
shapeConstr.bs <- function(x, shape, range = NULL, ...){

  # Extract specifications
  intercept <- attr(x, "intercept")
  ord <- attr(x, "degree") + 1
  bk <- attr(x, "Boundary.knots")
  ik <- attr(x, "knots")
  knots <- sort(c(rep(bk, ord), ik))

  # Create matrix
  mkShpCmat(shape, knots = knots, ord = ord, intercept = intercept,
    rng = range)
}
