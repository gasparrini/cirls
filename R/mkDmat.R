################################################################################
#
# Create a deriv matrix for shape-constrained splines with domain bounds
# Maybe export for "advanced" use?
#
################################################################################

mkDmat <- function(d, s, knots, ord, lower = -Inf, upper = Inf){

  #----- Create matrix

  # Create diff matrices for each derivative order
  alldm <- lapply(rev(seq_len(d)), dm, knots = knots, ord = ord)
  alldm <- c(alldm, list(diag(length(knots) - ord)))

  # Create full constraint matrix
  Dmat <- s * Reduce("%*%", alldm)

  #----- Remove rows with domain

  # Determine the constraints to be kept
  nk <- length(knots)
  upkn <- knots[(ord + 1):(nk - d)]
  lokn <- knots[(d + 1):(length(knots) - ord)]
  keep <- (upkn > lower) & (lokn < upper)

  # Remove from Cmat
  Dmat <- Dmat[keep,, drop = F]

  #----- Clean and return

  # Return
  Dmat
}
