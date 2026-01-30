test_that("crs basic estimation works", {
  set.seed(42)
  n <- 100
  x <- runif(n)
  z <- rbinom(n, 1, .5)
  y <- cos(2*pi*x) + z + rnorm(n, sd=.5)
  z <- factor(z)

  model <- crs(y ~ x + z, degree=3, segments=1, lambda=0.1, cv="none", kernel=TRUE, display.warnings=FALSE, display.nomad.progress=FALSE)
  
  expect_s3_class(model, "crs")
  expect_type(predict(model), "double")
  expect_equal(length(predict(model)), n)
})

test_that("crs works with default method", {
  set.seed(42)
  n <- 100
  x <- data.frame(x=runif(n))
  y <- cos(2*pi*x$x) + rnorm(n, sd=.5)
  
  model <- crs(x, y, degree=3, segments=1, cv="none", basis="additive", display.warnings=FALSE, display.nomad.progress=FALSE)
  expect_s3_class(model, "crs")
})
