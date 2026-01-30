test_that("clsd works", {
  set.seed(42)
  n <- 100
  x <- rnorm(n)
  
  # Compute density at a few points
  model <- clsd(x, degree=2, segments=1, xeval=seq(-2, 2, length=10), display.warnings=FALSE, display.nomad.progress=FALSE)
  
  expect_s3_class(model, "clsd")
  expect_type(model$density, "double")
  expect_length(model$density, 10)
  expect_true(all(model$density >= 0))
})
