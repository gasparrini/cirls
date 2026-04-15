################################################################################
#
# Method to create a constraint matrix for constraint on the bound of smooths
# strata method
#
################################################################################

#' @rdname boundConstr
#' @order 7
#' @export
boundConstr.strata <- function(x, ...){

  # Get info
  int <- attr(x, "intercept")
  ref <- attr(x, "ref")
  breaks <- attr(x, "breaks")
  ncat <- length(breaks) + 1

  # Get initial constraint matrix from default methods
  cm <- boundConstr.default(diag(ncat), intercept = TRUE, ...)

  # Adjust for reference: same as in dlnm:::strata
  if (!int & any(cm$Cmat[,1] != 0)) warning(
    "Need an intercept included to bound constrain on the left")
  if (ref > 0){
    cm$Cmat <- cm$Cmat[, -ref, drop = FALSE]
    if (int) cm$Cmat <- cbind(1, cm$Cmat)
  }

  # Return
  cm
}
