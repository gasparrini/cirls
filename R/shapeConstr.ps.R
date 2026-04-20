################################################################################
#
# Shape constraint matrix method:
# ps method
#
################################################################################

#' @rdname shapeConstr
#' @order 6
#' @export
shapeConstr.ps <- function(x, shape, ...){

  # Same as B-Splines
  shapeConstr.bs(x, shape, ...)
}
