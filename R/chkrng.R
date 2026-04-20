################################################################################
#
# Check range parameter for the shapeConstr method
#
################################################################################

chkrng <- function(rng, knots, msg = FALSE){

  # If not provided, use Inf
  if (is.null(rng)){
    rng <- c(-Inf, Inf)
  } else {
    # Otherwise, sort to avoid issues
    rng <- sort(rng)

    # Check the right length is provided
    # If only a single value, is interpreted as the max of the range
    if (length(rng) != 2){
      rng <- sort(c(rng, -Inf)[1:2])
      warning("'rng' should be a vector of length 2; using rng = ",
        deparse(rng))
    }

    # Check it contains any knot
    if (rng[1] > max(knots) | rng[2] < min(knots)){
      warning("No knot within provided 'rng'")
    } else {

      # If 'msg' is switched on, inform the user of range covered
      if (msg){
        ck <- c(max(knots[knots < rng[1]]), min(knots[knots > rng[2]]))
        message("Constraining range ", deparse(ck))
      }
    }
  }

  # Return
  rng
}
