################################################################################
#
# Method to create a constraint matrix for constraint on the bound of smooths
# bs method
#
################################################################################

# By default deg is set to the degree of the spline, for smoother convergence

#' @rdname boundConstr
#' @order 4
#' @export
boundConstr.bs <- function(x, deg = NULL, ...)
{

  # Get info from basis
  degree <- attr(x, "degree")
  df <- ncol(x)
  int <- attr(x, "intercept")

  # If necessary use degree of spline as default value
  deg <- deg %||% degree
  if (deg > degree) warning(paste0("Setting 'deg' greater than the degree",
    " of the spline will also constrain towards interior knots"))

  # Get constraint matrix
  boundConstr.default(x, deg = deg, intercept = int, ...)
}
