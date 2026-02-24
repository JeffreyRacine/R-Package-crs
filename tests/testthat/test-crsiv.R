test_that("crsiv works with Tikhonov", {
  set.seed(42)
  n <- 50 # smaller sample for speed
  v <- rnorm(n, sd=0.1)
  w_vec <- runif(n, -2, 2)
  y_vec <- w_vec + v + rnorm(n, sd=0.1)
  
  z <- data.frame(z=w_vec)
  w <- data.frame(w=w_vec)
  
  # crsiv returns a crsiv object (inheriting from crs)
  model <- crsiv(y=y_vec, z=z, w=w, method="Tikhonov", alpha=0.1, cv="none", basis="additive", display.warnings=FALSE, display.nomad.progress=FALSE)
  
  expect_s3_class(model, "crsiv")
  expect_s3_class(model, "crs")
})

test_that("crsiv works with Landweber-Fridman", {
  set.seed(42)
  n <- 50
  v <- rnorm(n, sd=0.1)
  w_vec <- runif(n, -2, 2)
  z_vec <- w_vec + v
  y_vec <- z_vec + v + rnorm(n, sd=0.1)
  
  z <- data.frame(z=z_vec)
  w <- data.frame(w=w_vec)
  
  # Use few iterations for speed
  model <- crsiv(y=y_vec, z=z, w=w, method="Landweber-Fridman", iterate.max=2, cv="none", basis="additive", display.warnings=FALSE, display.nomad.progress=FALSE)
  
  expect_s3_class(model, "crsiv")
  expect_s3_class(model, "crs")
  expect_type(model$phi, "double")
})

test_that("crsiv restores crs.messages option on internal error", {
  old_option <- getOption("crs.messages")
  on.exit(options(crs.messages = old_option), add = TRUE)

  options(crs.messages = TRUE)

  set.seed(42)
  n <- 25
  v <- rnorm(n, sd = 0.1)
  w_vec <- runif(n, -2, 2)
  z_vec <- w_vec + v
  y_vec <- z_vec + v + rnorm(n, sd = 0.1)

  z <- data.frame(z = z_vec)
  w <- data.frame(w = w_vec)

  expect_error(
    crsiv(
      y = y_vec,
      z = z,
      w = w,
      method = "Tikhonov",
      basis = "not-a-valid-basis",
      display.warnings = FALSE,
      display.nomad.progress = FALSE
    )
  )

  expect_identical(getOption("crs.messages"), TRUE)
})
