################################################################################
# Function to shape constrain cross-bases

#### Notes
# If a shapeConstr method, would allow to be called something like constr = ~ shape(cb)
# constr = ~ cb(cb) would be weirder I think
# If arguments left like that, would need to remove the argument "shape" from the generic method
# Plans to add the possibility to constrain at a given lag or a given at for lags
# Needs to work both for lags bs and ns

shapeConstr.crossbasis <- function(x, varshape, lagshape,
  overall = TRUE)
{

  #----- Initial checks

  # If any of the shape arguments are missing, use "pos"
  # Allows having the specificities of each basis
  if(missing(varshape)) varshape <- "pos"
  if(missing(lagshape)) lagshape <- "pos"

  #----- Recreate template bases for each dimensions
  # Used to rely on shapeConstr.onebasis

  # A simple basis for var dimension
  rng <- attr(x, "range")
  at <- seq(rng[1], rng[2], length.out = attr(x, "df")[1] + 2)
  varbasis <- do.call(dlnm::onebasis, c(list(x = at), attr(x, "argvar")))

  # Lag basis
  lags <- attr(x, "lag")
  lagvec <- seq(lags[1], lags[2], by = 1)
  lagbasis <- do.call(dlnm::onebasis, c(list(x = lagvec), attr(x, "arglag")))

  #----- Create constraints

  # In var dimension
  Cvar <- shapeConstr(varbasis, shape = varshape)

  # In lag dimension
  Clag <- shapeConstr(lagbasis, shape = lagshape)

  # Put together
  cmlist <- Map("%x%", Cvar, Clag)

  #----- Return
  cmlist
}
