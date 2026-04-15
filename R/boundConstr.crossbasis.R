################################################################################
#
# Method to create a constraint matrix for constraint on the bound of smooths
# crossbasis method
#
################################################################################

#' @export
boundConstr.crossbasis <- function(x, var = list(), lag = list()){

  # Extract info on crossbasis
  cbattr <- attributes(x)
  dims <- cbattr$df # Dimensions
  # Full ranges to prepare constraint matrix
  xrng <- cbattr$range
  lagrng <- cbattr$lag

  # Recreate template bases for each dimensions
  varbasis <- do.call(dlnm::onebasis, c(list(x = xrng), cbattr$argvar))
  lagbasis <- do.call(dlnm::onebasis, c(list(x = lagrng), cbattr$arglag))

  # Marginal constraint matrices
  Cvar <- do.call(boundConstr, c(list(x = varbasis), var))
  Clag <- do.call(boundConstr, c(list(x = lagbasis), var))

  # Put together
  Map("%x%", Cvar, Clag)
}
