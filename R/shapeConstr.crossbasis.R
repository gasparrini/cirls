################################################################################
# Function to shape constrain cross-bases

#' Constraints for distributed lag linear and nonlinear models
#'
#' @description
#' Methods to generate constraint matrices associated with [crossbasis][dlnm::crossbasis()] objects, allowing the fitting of constrained distributed-lag nonlinear models (DLNM). Designed for use within the [constr][buildCmat()] interface in `cirls`. **At the moment, constraining can only be performed on a single dimension at a time.**
#'
#' @param x A `crossbasis` object.
#' @param dim One of `"var"` or `"lag"`, the dimension on which the constraint will be applied. If `NULL`, the default is to apply the constraint on the `"var"` dimension if it is *nonlinear* (i.e. has more than one degree of freedom), or to apply it to the `"lag"` dimension otherwise (i.e. when a *linear* distributed lag model is specified).
#' @param slice A numeric vector of length 2 restricting the constraint on a specific range of the *other* dimension. By default no slicing is performed. See details.
#' @param overall Only when `dim = "var"`, logical indicating whether the constraint should be applied only on the overall cumulative association or across all lags. Note that `slice` can be used to define the constraint over a restricted cumulative lag range.
#' @param ... Parameters specific to the type of constraint to be applied. Includes for instance `shape` for [shapeConstr][shapeConstr()], or `value` for [boundConstr][boundConstr()]. See the main help page of the relevant generic.
#'
#' @details
#' Constraint matrices for `crossbasis` objects are built by generating the marginal constraint matrix for the requested dimension (through the `dim` argument) and expanding it through a Kronecker product with a (potentially modified) identity matrix of the dimension of the opposite dimension. This identity matrix can be modified to account for specificities of the basis for the other dimension (e.g. necessary for `ns`).
#'
#' ## Overall vs full surface
#'
#' By default, the constraint is applied over the whole surface. This means that if, say, the `var` dimension is shape-constrained, then the shape will hold true for any lag. When `overall` is switched to `TRUE` (only when `dim = "var"`), the constrain will hold for the *overall cumulative association* only (see [crossreduce][dlnm::crossreduce()]). This means the constraint could be violated at specific lags.
#'
#' ## Slicing
#'
#' The `slice` argument can be used to restrict the constraint on a specific range of the other dimension. For instance, when `dim = "var"`, setting `slice = c(0, 5)` (say) will enforce the constraint only over lags 0 to 5, meaning it could be violated at other lags (> 5 in this example). Note that the specific range for slicing depends on the specific basis for the opposite dimension. In the previous example, if the lag dimension is specified by splines, the constraint will be enforced between the closest knots that contain the range specified by `slice`. See the `range` argument in [shapeConstr][shapeConstr()] for additional details.
#'
#' When `overall = TRUE`, slicing will only constraint the cumulative association *over the range specific by `slice`*.
#'
#' ## Note
#'
#' These methods only allow to specify a single constraint on **one dimension at a time**. Specifying several constraints can be done by adding terms in the [constr][buildCmat()] formula.
#'
#' At the moment, **when the var dimension is nonlinear, it is not recommended to constrain the lag dimension**. This is because the centering performed by `crosspred` can impact the constraint definition.
#'
#' @returns A list containing the constraint matrix `Cmat`, and lower/upper bound vectors (`lb` and `ub`, respectively).
#'
#' @seealso [shapeConstr][shapeConstr()], [boundConstr][boundConstr()] for the generic method. [buildCmat][buildCmat()] for how to specify constraints.
#'
#' @references
#' Gasparrini, A., Armstrong, B., 2013. Reducing and meta-analysing estimates from distributed lag non-linear models. *BMC Medical Research Methodology* **13**, **1**. [DOI:10.1186/1471-2288-13-1](https://doi.org/10.1186/1471-2288-13-1)
#'
#' Gasparrini, A., Armstrong, B., Kenward, M.G., 2010. Distributed lag non-linear models. *Statistics in Medicine* **29**, **2224–2234**. [DOI:10.1002/sim.3940](https://doi.org/10.1002/sim.3940)
#'
#' @example inst/examples/ex_london_dlnm.R
#' @example inst/examples/ex_london_dlm.R
#'
#' @export
shapeConstr.crossbasis <- function(x, dim = NULL, slice = NULL, overall = FALSE,
  ...)
{
  # Call cbConstr
  cbConstr(x, constr = "shape", pars = list(...),
    dim = dim, slice = slice, overall = overall)
}
