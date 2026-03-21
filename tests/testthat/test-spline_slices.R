################################################################################
#
# Test shape constraints on subdomains
#
################################################################################

library(splines)

#----- Generate data ------

# Parameters
n <- 1000

# Generate a simple parabolic relationship
set.seed(1)
X <- seq(0, 1, length.out = n)
eta <- -.5 * sin(X * pi)
Y <- rnorm(n, eta, sd = .5)

#----- Spline basis -----

test_that("Slicing shape constraints work", {

  # Loop across bases
  for (bn in c("bs", "ns", "ps")){

    # Generate basis
    b <- do.call(bn, list(x = X, df = 10))

    # Fit and predict
    cm <- predict(glm(Y ~ b, method = "cirls.fit",
      constr = ~ shape(b, shape = "inc", range = c(.2, Inf))))

    # Test results
    expect_true(all(diff(cm[X > 0.2]) >= -sqrt(.Machine$double.eps)))
    expect_true(any(diff(cm[X < 0.1]) < -sqrt(.Machine$double.eps)))
  }

  # Now factor
  Xf <- cut(X, 10)
  cm <- predict(glm(Y ~ Xf, method = "cirls.fit",
    constr = ~ shape(Xf, shape = "ccv", range = c(3, Inf))))

  # Test results
  expect_true(all(diff(cm[X > 0.2]) >= -sqrt(.Machine$double.eps)))
  expect_true(any(diff(cm[X < 0.2]) < -sqrt(.Machine$double.eps)))
})

# # Plot
# plot(X, Y, col = "grey", pch = 16)
# lines(X, cm, col = 2, lwd = 3)
# abline(v = knots, lty = 3)

