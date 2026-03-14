################################################################################
#
# Tests for shape-constrained DLNM
#
################################################################################

library(dlnm)

#-----------------------------
# Simulate data
#-----------------------------

# Some parameters
maxlag <- 21

#----- Functions for simulation

# On variable dimension
fpara <- function(x) -.5 * sin(x * pi)

# Lag dimension
wdecay <- function(lag) exp(-lag/5)

# Full surface
fsurf <- function(x,lag) 0.1 * fpara(x) * wdecay(lag)

#----- Simulate data

# Use actual data for X
X <- london$tmean
X <- (X - min(X, na.rm = T)) / diff(range(X, na.rm = T))

# Compute cumulated effect
seqlag <- 0:maxlag
Q <- sapply(seqlag, function(l) c(rep(NA, l), X[1:(length(X) - l)]))
cumeff <- apply(Q, 1, function(hist) sum(fsurf(hist, seqlag)))

# Simulate response
set.seed(21)
suppressWarnings(y <- rpois(length(X), exp((log(20) + cumeff))))

#-----------------------------
# test shape-constrained DLNM
#-----------------------------

test_that("monotone-constr dlnm works", {

# List of type of basis functions
funlist <- c("strata", "bs", "ns")

# Loop over the types
for (vf in funlist) for (lf in funlist) for (shp in c("inc", "dec")){
  # Define crossbasis
  cb <- crossbasis(X, lag = maxlag, argvar = list(fun = vf, df = 10),
    arglag = list(fun = lf, df = 5))

  # Constraint matrix and fit model
  Cmat <- shapeConstr.crossbasis(cb, varshape = shp)
  cmodel <- glm(y ~ cb, family = "quasipoisson",
    method = "cirls.fit", Cmat = list(cb = Cmat$Cmat))

  # Prediction - each column (lag) should be non-decreasing
  ccp <- crosspred(cb, cmodel)
  if (shp == "inc") expect_true(
    all(diff(ccp$matfit) >= -sqrt(.Machine$double.eps)))
  if (shp == "dec") expect_true(
    all(diff(ccp$matfit) <= sqrt(.Machine$double.eps)))
}

})


#-----------------------------
# test shape-constrained cumulative overall
#-----------------------------

cb <- crossbasis(X, lag = maxlag, argvar = list(fun = "bs", df = 10),
    arglag = list(fun = "ns", df = 5))
umodel <- glm(y ~ cb, family = "quasipoisson")
cp <- crosspred(cb, umodel)
plot(cp, ptype = "overall")
