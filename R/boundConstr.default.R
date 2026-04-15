################################################################################
#
# Bound constraint
# default method
#
################################################################################

#' @rdname boundConstr
#' @order 2
#' @export
boundConstr.default <- function(x, value = 0, side = c("right", "left", "both"),
  deg = NULL, intercept = FALSE, ...)
{

  # Extract info
  df <- ncol(x)

  # Check side parameter
  side <- match.arg(side)

  # Check degree parameters
  deg <- deg %||% 1
  if (deg > df | deg < 1) stop(
      "'deg' should be an integer between 1 and the number of bases")

  # Check intercept
  if (!intercept & side != "right") warning(
    "Need an intercept included to bound constrain on the left")

  # Create Cmat depending on side
  cr <- cl <- NULL
  if (side != "right"){
    degl <- deg - 1 + intercept
    cl <- cbind(diag(degl), matrix(0, degl, df - degl))
  }
  if (side != "left"){
    cr <- cbind(matrix(0, deg, df - deg), diag(deg))
  }
  Cmat <- rbind(cr, cl)

  # Return with bounds
  list(Cmat = Cmat, lb = rep(value, NROW(Cmat)), ub = rep(value, NROW(Cmat)))
}
