###############################################################################
# Non-negative coefficient example with the london dataset

library(splines)

### Association between Ozone and mortality
# Nonnegative constraint on Ozone
model <- glm(death ~ o3, data = london, family = "quasipoisson",
  method = "cirls.fit", cons = ~ shape(o3, "pos"))

# Coefficient and confidence interval
coef(model)[2]
confint(model)[2,]

# Comparing to an unconstrained model
umodel <- uncons(model)
coef(umodel)[2]
confint(umodel)[2,]
