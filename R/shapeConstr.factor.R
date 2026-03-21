################################################################################
#
# Shape constraint matrix method:
# Default method
#
################################################################################

#' @rdname shapeConstr
#' @export
shapeConstr.factor <- function(x, shape, range = NULL, intercept = FALSE, ...) {

  # Get the design matrix
  xmat <- model.matrix(~ x)

  # Get initial constraint matrix from default methods
  Cmat <- shapeConstr.default(xmat, shape = shape, range = range,
    intercept = TRUE)$Cmat

  # Apply contrast if intercept is not included
  # NB: no contrast is applied in R when the model does not include an intercept
  if (isFALSE(intercept)){
    ctr <- stats::contrasts(x)
    Cmat <- Cmat %*% ctr
  }

  # Bounds
  lb <- rep(0, NROW(Cmat))
  ub <- rep(Inf, NROW(Cmat))

  # Return
  list(Cmat = Cmat, lb = lb, ub = ub)
}
