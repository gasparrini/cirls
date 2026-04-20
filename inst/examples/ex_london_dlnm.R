###############################################################################
# Constrained DLNM: var dimension

library(dlnm)
library(splines)

#----------------
# Association between temperature and mortality
#----------------

# Number of years and dow
ny <- length(unique(format(london$date, "%Y")))

# Create a flexible crossbasis
cb <- crossbasis(london$tmean, lag = 21,
  argvar = list(fun = "bs", degree = 2, df = 10),
  arglag = list(fun = "ns", knots = logknots(21, df = 5)))

#----- Unconstrained model

# Fit the model
um <- glm(age0_64 ~ cb + ns(date, df = 7 * ny) + dow, data = london,
  family = "quasipoisson")
ucp <- crosspred(cb, um, cen = 20)

# Plot overall
plot(ucp, ptype = "overall", lwd = 2)

# Plot slices
plot(ucp, ptype = "slice", lag = 0, lwd = 2)
lines(ucp, ptype = "slice", lag = 5, col = 2, lwd = 2)
lines(ucp, ptype = "slice", lag = 15, col = 3, lwd = 2)
legend("topleft", legend = sprintf("lag = %i", c(0, 5, 15)), col = 1:3,
  lty = 1)

#----- Constrained DLNM

# Fit the model with convexity constraint
scm <- glm(age0_64 ~ cb + ns(date, df = 7 * ny) + dow, data = london,
  family = "quasipoisson", method = "cirls.fit",
  constr = ~ shape(cb, shape = "cvx"))
scp <- crosspred(cb, scm, cen = 20)

# Plot slices
plot(scp, ptype = "slice", lag = 0, lwd = 2)
lines(scp, ptype = "slice", lag = 5, col = 2, lwd = 2)
lines(scp, ptype = "slice", lag = 15, col = 3, lwd = 2)
legend("topleft", legend = sprintf("lag = %i", c(0, 5, 15)), col = 1:3,
  lty = 1)

#----- Overall only

# Fit the model with a constraint on overall only
ocm <- glm(age0_64 ~ cb + ns(date, df = 7 * ny) + dow, data = london,
  family = "quasipoisson", method = "cirls.fit",
  constr = ~ shape(cb, shape = "cvx", overall = TRUE))
ocp <- crosspred(cb, ocm, cen = 20)

# Plot slices: they don't necessarily respect the constraints
plot(ocp, ptype = "slice", lag = 0, lwd = 2)
lines(ocp, ptype = "slice", lag = 5, col = 2, lwd = 2)
lines(ocp, ptype = "slice", lag = 15, col = 3, lwd = 2)
legend("topleft", legend = sprintf("lag = %i", c(0, 5, 15)), col = 1:3,
  lty = 1)

# Plot overall: constrained
plot(ocp, ptype = "overall", lwd = 2)
lines(ucp, ptype = "overall", col = 2, lwd = 2, ci = "lines")
legend("topleft", legend = c("Constrained", "Unconstrained"),
  col = 1:2, lty = 1:2)

#----- Slices

# Fit model
subcm <- glm(age0_64 ~ cb + ns(date, df = 7 * ny) + dow, data = london,
  family = "quasipoisson", method = "cirls.fit",
  constr = ~ shape(cb, shape = "cvx", slice = c(0, 1)))
subcp <- crosspred(cb, subcm, cen = 20)

# Plot slices: lags 0 and 5 are constrained but not lag 15
plot(subcp, ptype = "slice", lag = 0, lwd = 2)
lines(subcp, ptype = "slice", lag = 5, col = 2, lwd = 2)
lines(subcp, ptype = "slice", lag = 15, col = 3, lwd = 2)
legend("topleft", legend = sprintf("lag = %i", c(0, 5, 15)), col = 1:3, lty = 1)


