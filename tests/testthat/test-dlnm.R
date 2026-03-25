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
wpeak <- function(lag) 5 * dnorm(lag, 5, 5)

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
for (vf in funlist) for (lf in funlist){
  # Define crossbasis
  cb <- crossbasis(X, lag = maxlag, argvar = list(fun = vf, df = 10),
    arglag = list(fun = lf, df = 5))

  # Constraint matrix and fit model
  cmodel <- glm(y ~ cb, family = "quasipoisson",
    method = "cirls.fit", constr = ~ shape(cb, vshape = "inc"))

  # Prediction - each column (lag) should be non-decreasing
  ccp <- crosspred(cb, cmodel, cen = .5)
  expect_true(all(diff(ccp$matfit) >= -sqrt(.Machine$double.eps)))
}

})


#-----------------------------
# test shape-constrained cumulative overall
#-----------------------------


test_that("overall constraining works", {

# List of type of basis functions
funlist <- c("strata", "bs", "ns")

# Loop over the types
for (vf in funlist) for (lf in funlist){
  # Define crossbasis
  cb <- crossbasis(X, lag = maxlag, argvar = list(fun = vf, df = 10),
    arglag = list(fun = lf, df = 5))

  # Constraint matrix and fit model
  cmodel <- glm(y ~ cb, family = "quasipoisson",
    method = "cirls.fit", constr = ~ shape(cb, vshape = "inc", overall = TRUE))

  # Prediction - each column (lag) should be non-decreasing
  ccp <- crosspred(cb, cmodel, cen = .5)
  expect_true(all(diff(ccp$allfit) >= -sqrt(.Machine$double.eps)))

  # Plot
  # plot(ccp, ptype = "overall")
  # plot(ccp, ptype = "slices", lag = 0)
}

})

#-----------------------------
# test slices
#-----------------------------

test_that("slicing works", {

# List of type of basis functions
funlist <- c("strata", "bs", "ns")

# Loop over the types
for (vf in funlist) for (lf in funlist){
  # Define crossbasis
  cb <- crossbasis(X, lag = maxlag, argvar = list(fun = vf, df = 10),
    arglag = list(fun = lf, df = 5))

  # Constraint matrix and fit model
  cmodel <- glm(y ~ cb, family = "quasipoisson", method = "cirls.fit",
    constr = ~ shape(cb, vshape = "inc", lrange = c(0, 5)))

  # Prediction - each column (lag) should be non-decreasing
  ccp <- crosspred(cb, cmodel, cen = .5)
  expect_true(all(diff(ccp$matfit[,1]) >= -sqrt(.Machine$double.eps)))
  expect_true(all(diff(ccp$matfit[,5]) >= -sqrt(.Machine$double.eps)))
  # expect_false(all(diff(ccp$matfit[,15]) >= -sqrt(.Machine$double.eps)))

  # Plot
  # plot(ccp, ptype = "overall")
  # plot(ccp, ptype = "slices", lag = 0)
}

})

################################################################################
# DLM test

#--------------------------
# Resimulate data
#--------------------------

# Create cumulative effect
cumeff2 <- Q %*% (wpeak(0:maxlag) / 5)

# Simulate response
set.seed(24)
suppressWarnings(y <- rpois(length(X), exp((log(20) + cumeff2))))

#----- Unconstrained DLM

# Basis
dlb <- crossbasis(X, lag = maxlag, argvar = list(fun = "lin"),
    arglag = list(fun = "strata", df = 5))

# Fit model
udlm <- glm(y ~ dlb, family = "quasipoisson")

# Look at it
crosspred(dlb, udlm, cen = .5) |> plot(ptype = "slices", var = 1)


#----- Constrained DLM

# Fit model
cdlm <- glm(y ~ dlb, family = "quasipoisson",
    method = "cirls.fit", constr = ~ shape(dlb, lshape = "dec"))

# Look at it
crosspred(dlb, cdlm, cen = .5) |> plot(ptype = "slices", var = 1)
