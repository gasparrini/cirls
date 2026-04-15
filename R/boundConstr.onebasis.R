################################################################################
#
# Method to create a constraint matrix for constraint on the bound of smooths
# onebasis method
#
################################################################################

#' @rdname boundConstr
#' @order 8
#' @export
boundConstr.onebasis <- function(x, ...){

  # Extract the right method
  fun <- attr(x, "fun")
  met <- paste0("boundConstr.", fun)
  if (!met %in% utils::methods("boundConstr")) {
    warning(paste0("No existing 'boundConstr' method for '", fun,
      "' functions. Using default method."))
    met <- "boundConstr.default"
  }

  # Call the right method
  pars <- utils::modifyList(list(x = x), list(...))
  do.call(met, pars)
}
