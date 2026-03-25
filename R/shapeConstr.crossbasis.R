################################################################################
# Function to shape constrain cross-bases

#### Notes
# If a shapeConstr method, would allow to be called something like constr = ~ shape(cb)
# constr = ~ cb(cb) would be weirder I think
# If arguments left like that, would need to remove the argument "shape" from the generic method
# Plans to add the possibility to constrain at a given lag or a given at for lags
# Needs to work both for lags bs and ns

#' @exportS3Method
shapeConstr.crossbasis <- function(x, vshape, lshape,
  overall = FALSE, vrange = NULL, lrange = NULL, ...)
{
  # Extract info on crossbasis
  cbattr <- attributes(x)
  dims <- cbattr$df # Dimensions
  # Full ranges to prepare constraint matrix
  xrng <- cbattr$range
  lag <- cbattr$lag

  #----- Initial checks

  # If any of the shape arguments are missing, use "pos"
  # Allows having the specificity of each basis
  if(missing(vshape)) vshape <- "pos"
  if(missing(lshape)) lshape <- "pos"

  # Check range and lag to be constrained
  vrange <- chkrng(vrange, xrng)
  lrange <- chkrng(lrange, lag)

  #----- Recreate template bases for each dimensions
  # Used to rely on shapeConstr.onebasis

  # A simple basis for var dimension
  at <- seq(xrng[1], xrng[2], length.out = dims[1] + 2)
  varbasis <- do.call(dlnm::onebasis, c(list(x = at), cbattr$argvar))

  # Lag basis
  lagvec <- seq(lag[1], lag[2], by = 1)
  lagbasis <- do.call(dlnm::onebasis, c(list(x = lagvec), cbattr$arglag))

  #----- Create marginal constraint matrices

  # In var dimension
  Cvar <- shapeConstr(varbasis, shape = vshape, range = vrange)

  # In lag dimension
  Clag <- shapeConstr(lagbasis, shape = lshape, range = lrange)

  #----- Final constraint matrix

  # Constrain overall only
  if (isTRUE(overall)){

    # This is the same M matrix as in Gasparrini Armstrong 2013
    M <- diag(attr(x, "df")[1]) %x% t(colSums(lagbasis))
    ovCmat <- Cvar$Cmat %*% M
    cmlist <- utils::modifyList(Cvar, list(Cmat = ovCmat))

  } else {
    # Or the whole surface: a simple kronecker product
    cmlist <- Map("%x%", Cvar, Clag)
  }

  #----- Return
  cmlist
}
