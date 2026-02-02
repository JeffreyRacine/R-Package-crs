context("crsivderiv")

test_that("crsivderiv works with Landweber-Fridman", {
  set.seed(42)
  n <- 100
  v <- rnorm(n, sd=0.1)
  w_vec <- runif(n, -2, 2)
  z_vec <- 0.2*w_vec + v
  y_vec <- z_vec^2 - 0.5*v + rnorm(n, sd=0.05)

  z <- data.frame(z=z_vec)
  w <- data.frame(w=w_vec)

  # Basic execution test
  # Use small iterate.max for speed in tests
  model <- crsivderiv(y=y_vec, z=z, w=w, iterate.max=5, cv="none", basis="additive", display.nomad.progress=FALSE)
  
  expect_s3_class(model, "crsivderiv")
  expect_s3_class(model, "crs")
  expect_type(model$phi, "double")
  expect_type(model$phi.prime, "double")
  expect_equal(length(model$phi), n)
  expect_equal(length(model$phi.prime), n)
  expect_equal(ncol(model$phi.mat), model$num.iterations)
})

test_that("crsivderiv handles evaluation data correctly", {
  set.seed(42)
  n <- 50
  w_vec <- rnorm(n)
  z_vec <- 0.2*w_vec + rnorm(n, sd=0.1)
  y_vec <- z_vec^2 + rnorm(n, sd=0.1)
  
  z <- data.frame(z=z_vec)
  w <- data.frame(w=w_vec)
  
  zeval <- data.frame(z=seq(min(z_vec), max(z_vec), length.out=n))
  weval <- data.frame(w=rep(0, n))
  
  model <- crsivderiv(y=y_vec, z=z, w=w, zeval=zeval, weval=weval, iterate.max=2, cv="none", basis="additive", display.nomad.progress=FALSE)
  
  expect_equal(length(model$phi), n)
  expect_equal(length(model$phi.prime), n)
})

test_that("crsivderiv works with exogenous predictors", {
  set.seed(42)
  n <- 50
  x_vec <- rnorm(n)
  w_vec <- rnorm(n)
  z_vec <- 0.2*w_vec + 0.1*x_vec + rnorm(n, sd=0.1)
  y_vec <- z_vec^2 + 0.5*x_vec + rnorm(n, sd=0.1)
  
  z <- data.frame(z=z_vec)
  w <- data.frame(w=w_vec)
  x <- data.frame(x=x_vec)
  
  model <- crsivderiv(y=y_vec, z=z, w=w, x=x, iterate.max=2, cv="none", basis="additive", display.nomad.progress=FALSE)
  
  expect_equal(length(model$phi), n)
  expect_type(model$phi, "double")
})
