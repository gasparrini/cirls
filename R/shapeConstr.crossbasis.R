################################################################################
# Function to shape constrain cross-bases

#### Notes
# If a shapeConstr method, would allow to be called something like constr = ~ shape(cb)
# constr = ~ cb(cb) would be weirder I think
# If arguments left like that, would need to remove the argument "shape" from the generic method
# Plans to add the possibility to constrain at a given lag or a given at for lags
# Needs to work both for lags bs and ns

#' @exportS3Method
shapeConstr.crossbasis <- function(x, dim = NULL, slice = NULL, overall = FALSE,
  ...)
{
  # Call cbConstr
  cbConstr(x, constr = "shape", pars = list(...),
    dim = dim, slice = slice, overall = overall)
}
