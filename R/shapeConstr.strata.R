################################################################################
#
# Shape constraint matrix method:
# Strata
#
################################################################################

#' @rdname shapeConstr
#' @order 8
#' @export
shapeConstr.strata <- function(x, shape, range = NULL, ...){

  # Get parameters
  ref <- attr(x, "ref")
  intercept <- attr(x, "intercept")
  breaks <- attr(x, "breaks")
  ncat <- length(breaks) + 1

  # Define range in terms of variables
  rng <- chkrng(range, c(-Inf, breaks, Inf), msg = FALSE)
  keep <- (c(-Inf, breaks) < rng[2]) & (c(breaks, Inf) > rng[1])
  rng <- range(which(keep))

  # Call the default method to get first matrix (always with intercept)
  cmlist <- shapeConstr.default(diag(ncat), shape = shape, range = rng,
    intercept = TRUE, ...)

  # Adjust for reference
  Cmat <- cmlist$Cmat
  if (ref > 0){
    Cmat <- Cmat[, -ref, drop = FALSE]
    if (intercept){
      cp <- chkshp(shape)
      Cmat <- if (cp[[1]][1] == 0) cbind(1, Cmat) else cbind(0, Cmat)
    }
  }

  # Recompute bounds
  lb <- rep(0, NROW(Cmat))
  ub <- rep(Inf, NROW(Cmat))

  # Return
  list(Cmat = Cmat, lb = lb, ub = ub)
}
