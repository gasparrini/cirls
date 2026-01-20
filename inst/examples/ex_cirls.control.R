###############################################################################
# Examples of control parameters

# The simplest way to specify parameters is through ...
model <- glm(death ~ o3, data = london, family = "quasipoisson",
  method = "cirls.fit", cons = ~ shape(o3, "pos"))
osqpmodel <- glm(death ~ o3, data = london, family = "quasipoisson",
  method = "cirls.fit", cons = ~ shape(o3, "pos"), qp_solver = "osqp")

# Comparing the control
model$control
osqpmodel$control

# Alternatively through the control argument
osqpmodel2 <- glm(death ~ o3, data = london, family = "quasipoisson",
  method = "cirls.fit", control = list(
    cons = ~ shape(o3, "pos"), qp_solver = "osqp"))
all.equal(osqpmodel$control, osqpmodel2$control)

# However, both cannot be used at the same time
osqpmodel3 <- glm(death ~ o3, data = london, family = "quasipoisson",
  method = "cirls.fit", control = list(cons = ~ shape(o3, "pos")),
  qp_solver = "osqp")
all.equal(osqpmodel$control, osqpmodel3$control)

# The control argument takes precedence
model2 <- glm(death ~ o3, data = london, family = "quasipoisson",
  method = "cirls.fit", control = list(cons = ~ shape(o3, "pos")),
  cons = ~ shape(o3, "neg"))
model2$control$cons
