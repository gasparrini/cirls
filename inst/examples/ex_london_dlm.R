###############################################################################
# Constrained distributed lag models (DLM)

library(dlnm)
library(splines)

#----------------
# Association between PM2.5 and mortality
#----------------

# Number of years and dow
ny <- length(unique(format(london$date, "%Y")))

# Create a flexible crossbasis
cb <- crossbasis(london$pm10, lag = 15, argvar = list(fun = "lin"),
  arglag = list(fun = "ns", knots = logknots(15, df = 4)))

#----- Unconstrained model

# Fit the model
um <- glm(age0_64 ~ cb + ns(date, df = 7 * ny) + dow, data = london,
  family = "quasipoisson")

# Plot lag-response function
ucp <- crosspred(cb, um, cen = 0, at = 10)
plot(ucp, ptype = "slices", lwd = 2, var = 10)

#----- Constrained to converge to zero

# Fit with bound constraint
# By default dim = "lag" (because var is linear) and value = 0
bm <- glm(age0_64 ~ cb + ns(date, df = 7 * ny) + dow, data = london,
  family = "quasipoisson", method = "cirls.fit", constr = ~ bound(cb))

# Plot
bcp <- crosspred(cb, bm, cen = 0, at = 10)
plot(bcp, ptype = "slices", lwd = 2, var = 10)

#----- Constrained to always be non-negative

# Fit with bound constraint
# By default dim = "lag" (because var is linear)
nnm <- glm(age0_64 ~ cb + ns(date, df = 7 * ny) + dow, data = london,
  family = "quasipoisson", method = "cirls.fit",
  constr = ~ shape(cb, shape = "pos"))

# Plot
nncp <- crosspred(cb, nnm, cen = 0, at = 10)
plot(nncp, ptype = "slices", lwd = 2, var = 10)

#----- Constrained to always decrease

# Fit with bound constraint
# By default dim = "lag" (because var is linear)
dm <- glm(age0_64 ~ cb + ns(date, df = 7 * ny) + dow, data = london,
  family = "quasipoisson", method = "cirls.fit",
  constr = ~ shape(cb, shape = "dec"))

# Plot
dcp <- crosspred(cb, dm, cen = 0, at = 10)
plot(dcp, ptype = "slices", lwd = 2, var = 10)

#----- Both decrease and towards zero

# Fit with bound constraint
# By default dim = "lag" (because var is linear)
dbm <- glm(age0_64 ~ cb + ns(date, df = 7 * ny) + dow, data = london,
  family = "quasipoisson", method = "cirls.fit",
  constr = ~ shape(cb, shape = "dec") + bound(cb))

# Plot
dbcp <- crosspred(cb, dbm, cen = 0, at = 10)
plot(dbcp, ptype = "slices", lwd = 2, var = 10)
