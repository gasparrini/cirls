################################################################################
#
# Shape constraint matrix method:
# linear term
# For convenience with the `dlnm` package
#
################################################################################

#' @rdname shapeConstr
#' @order 9
#' @export
shapeConstr.lin <- function(x, shape, ...){

  # Check shape parameter
  cpars <- chkshp(shape, 2)
  if (length(cpars) > 1) stop("Inconsistent shapes provided. ",
    "Linear terms only allow one of 'pos'/'inc' or 'neg'/dec'.")

  # Return only sign
  list(Cmat = matrix(cpars[[1]][2], 1, 1), lb = 0, ub = Inf)
}
