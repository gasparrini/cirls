###############################################################################
# Bound constraints example: cyclic seasonality of mortality

library(splines)

# Spline basis for seasonality
sbasis <- bs(london$doy, df = 7)

# Estimation of trend and cyclic seasonality with boundaries constrained to 0
model <- glm(death ~ date + sbasis, data = london, family = "quasipoisson",
  method = "cirls.fit", cons = ~ bound(sbasis))

# Reconstruct seasonality
sb <- do.call(bs, c(list(x = 1:366),
  attributes(sbasis)[c("degree", "knots", "Boundary.knots")]))
seas <- sb %*% coef(model)[-(1:2)]

# Plot
plot(seas, type = "l", xlab = "Day of year", ylab = "Seasonnality",
  ylim = c(-0.5, 0.1))

# With `bs`, deg can be reduced to obtain a less smooth convergence
model2 <- glm(death ~ date + sbasis, data = london, family = "quasipoisson",
  method = "cirls.fit", cons = ~ bound(sbasis, deg = 1))
seas2 <- sb %*% coef(model2)[-(1:2)]
lines(seas2, col = 2)
