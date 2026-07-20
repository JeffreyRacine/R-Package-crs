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
  model <- crsivderiv(y=y_vec, z=z, w=w, iterate.max=5, cv="none", basis="additive", display.nomad.progress=FALSE, display.warnings=FALSE)
  
  expect_s3_class(model, "crsivderiv")
  expect_s3_class(model, "crs")
  expect_type(model$phi, "double")
  expect_type(model$phi.prime, "double")
  expect_equal(length(model$phi), n)
  expect_equal(length(model$phi.prime), n)
  expect_equal(ncol(model$phi.mat), length(model$norm.stop))
  expect_equal(ncol(model$phi.prime.mat), length(model$norm.stop))
  expect_true(model$num.iterations >= 1L)
  expect_true(model$num.iterations <= length(model$norm.stop))
  expect_equal(model$phi, model$phi.mat[, model$num.iterations])
  expect_equal(model$phi.prime,
               model$phi.prime.mat[, model$num.iterations])
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
  
  model <- crsivderiv(y=y_vec, z=z, w=w, zeval=zeval, weval=weval, iterate.max=2, cv="none", basis="additive", display.nomad.progress=FALSE, display.warnings=FALSE)
  
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
  
  model <- crsivderiv(y=y_vec, z=z, w=w, x=x, iterate.max=2, cv="none", basis="additive", display.nomad.progress=FALSE, display.warnings=FALSE)
  
  expect_equal(length(model$phi), n)
  expect_type(model$phi, "double")
})

test_that("crsivderiv restores crs.messages option on internal error", {
  old_option <- getOption("crs.messages")
  on.exit(options(crs.messages = old_option), add = TRUE)

  options(crs.messages = TRUE)

  set.seed(42)
  n <- 30
  w_vec <- rnorm(n)
  z_vec <- 0.2 * w_vec + rnorm(n, sd = 0.1)
  y_vec <- z_vec^2 + rnorm(n, sd = 0.1)

  z <- data.frame(z = z_vec)
  w <- data.frame(w = w_vec)

  expect_error(
    crsivderiv(
      y = y_vec,
      z = z,
      w = w,
      basis = "not-a-valid-basis",
      iterate.max = 2,
      display.warnings = FALSE,
      display.nomad.progress = FALSE
    )
  )

  expect_identical(getOption("crs.messages"), TRUE)
})

test_that("crsivderiv works when crs.messages option is unset", {
  old_option <- getOption("crs.messages")
  on.exit(options(crs.messages = old_option), add = TRUE)

  options(crs.messages = NULL)

  set.seed(42)
  n <- 40
  w_vec <- rnorm(n)
  z_vec <- 0.2 * w_vec + rnorm(n, sd = 0.1)
  y_vec <- z_vec^2 + rnorm(n, sd = 0.1)

  z <- data.frame(z = z_vec)
  w <- data.frame(w = w_vec)

  model <- crsivderiv(
    y = y_vec,
    z = z,
    w = w,
    iterate.max = 2,
    cv = "none",
    basis = "additive",
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )

  expect_s3_class(model, "crsivderiv")
  expect_null(getOption("crs.messages"))
})

test_that("crsivderiv NOMAD progress keeps ownership at the IV wrapper", {
  test_time_values <- function(values) {
    i <- 0L
    force(values)
    function() {
      i <<- min(i + 1L, length(values))
      values[[i]]
    }
  }

  set.seed(23)
  n <- 30L
  v <- rnorm(n, sd = 0.1)
  w_vec <- runif(n, -2, 2)
  z_vec <- 0.2 * w_vec + v
  y_vec <- z_vec^2 - 0.5 * v + rnorm(n, sd = 0.05)
  z <- data.frame(z = z_vec)
  w <- data.frame(w = w_vec)

  traced <- capture_crs_progress_shadow_trace(
    suppressWarnings(
      crsivderiv(
        y = y_vec,
        z = z,
        w = w,
        iterate.max = 2,
        cv = "nomad",
        cv.threshold = 0,
        basis = "additive",
        degree.max = 2,
        segments.max = 2,
        nmulti = 2,
        max.bb.eval = 20,
        display.warnings = FALSE,
        display.nomad.progress = TRUE
      )
    ),
    force_renderer = "legacy",
    now = test_time_values(seq(0, 40, by = 0.25)),
    interactive = TRUE
  )

  lines <- unique(vapply(traced$trace, `[[`, character(1L), "line"))
  expect_true(any(grepl("Computing optimal smoothing", lines, fixed = TRUE)))
  expect_false(any(grepl("^\\[crs\\] Selecting spline model\\.\\.\\.", lines)))
  expect_false(any(grepl("multistart 1/2", lines, fixed = TRUE)))
  expect_false(any(grepl("fv=", lines, fixed = TRUE)))
  expect_false(any(grepl("Calling NOMAD", lines, fixed = TRUE)))
})

test_that("crsivderiv NOMAD progress does not change fixed-seed results", {
  set.seed(29)
  n <- 30L
  v <- rnorm(n, sd = 0.1)
  w_vec <- runif(n, -2, 2)
  z_vec <- 0.2 * w_vec + v
  y_vec <- z_vec^2 - 0.5 * v + rnorm(n, sd = 0.05)
  z <- data.frame(z = z_vec)
  w <- data.frame(w = w_vec)

  silent <- suppressWarnings(
    crsivderiv(
      y = y_vec,
      z = z,
      w = w,
      iterate.max = 2,
      cv = "nomad",
      cv.threshold = 0,
      basis = "additive",
      degree.max = 2,
      segments.max = 2,
      nmulti = 2,
      max.bb.eval = 20,
      random.seed = 41,
      display.warnings = FALSE,
      display.nomad.progress = FALSE
    )
  )

  noisy <- suppressWarnings(
    crsivderiv(
      y = y_vec,
      z = z,
      w = w,
      iterate.max = 2,
      cv = "nomad",
      cv.threshold = 0,
      basis = "additive",
      degree.max = 2,
      segments.max = 2,
      nmulti = 2,
      max.bb.eval = 20,
      random.seed = 41,
      display.warnings = FALSE,
      display.nomad.progress = TRUE
    )
  )

  expect_equal(as.numeric(noisy$phi), as.numeric(silent$phi), tolerance = 1e-12)
  expect_equal(as.numeric(noisy$phi.prime), as.numeric(silent$phi.prime), tolerance = 1e-12)
  expect_identical(noisy$num.iterations, silent$num.iterations)
})
