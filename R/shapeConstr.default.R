################################################################################
#
# Shape constraint matrix method:
# Default method
#
################################################################################

#' @rdname shapeConstr
#' @export
shapeConstr.default <- function(x, shape, intercept = FALSE, ...) {

  # Matrix dimension
  if (length(dim(x)) < 2) x <- as.matrix(x)
  ord <- ncol(x)

  # Check parameters
  cpars <- chkshp(shape, ord)

  # Create constraint matrices
  knots <- seq_len(2 * ord + !intercept)
  Cmat <- lapply(cpars, function(cp) dmat(cp[1], cp[2], knots, ord))
  Cmat <- do.call(rbind, Cmat)

  # Put together and remove redundant constraints
  if (!intercept) Cmat <- Cmat[,-1]
  chkc <- checkCmat(Cmat)
  if (length(chkc$redundant) > 0) Cmat <- Cmat[-chkc$redundant, , drop = F]

  # Return
  Cmat
}
