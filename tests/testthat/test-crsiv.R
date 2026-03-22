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

test_that("crsiv works when crs.messages option is unset", {
  old_option <- getOption("crs.messages")
  on.exit(options(crs.messages = old_option), add = TRUE)

  options(crs.messages = NULL)

  set.seed(42)
  n <- 40
  v <- rnorm(n, sd = 0.1)
  w_vec <- runif(n, -2, 2)
  z_vec <- w_vec + v
  y_vec <- z_vec + v + rnorm(n, sd = 0.1)

  z <- data.frame(z = z_vec)
  w <- data.frame(w = w_vec)

  model <- crsiv(
    y = y_vec,
    z = z,
    w = w,
    method = "Landweber-Fridman",
    iterate.max = 2,
    cv = "none",
    basis = "additive",
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )

  expect_s3_class(model, "crsiv")
  expect_null(getOption("crs.messages"))
})

test_that("crsiv NOMAD progress hands off to nested crs selectors", {
  test_time_values <- function(values) {
    i <- 0L
    force(values)
    function() {
      i <<- min(i + 1L, length(values))
      values[[i]]
    }
  }

  set.seed(17)
  n <- 30L
  v <- rnorm(n, sd = 0.1)
  w_vec <- runif(n, -2, 2)
  z_vec <- w_vec + v
  y_vec <- z_vec + v + rnorm(n, sd = 0.1)
  z <- data.frame(z = z_vec)
  w <- data.frame(w = w_vec)

  traced <- capture_crs_progress_shadow_trace(
    suppressWarnings(
      crsiv(
        y = y_vec,
        z = z,
        w = w,
        method = "Landweber-Fridman",
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
  expect_true(any(grepl("IV regression", lines, fixed = TRUE)))
  expect_true(any(grepl("^\\[crs\\] Selecting spline model\\.\\.\\.", lines)))
  expect_true(any(grepl("multistart 1/2", lines, fixed = TRUE)))
  expect_true(any(grepl("fv=", lines, fixed = TRUE)))
  expect_false(any(grepl("Calling NOMAD", lines, fixed = TRUE)))
})

test_that("crsiv NOMAD progress does not change fixed-seed results", {
  set.seed(19)
  n <- 30L
  v <- rnorm(n, sd = 0.1)
  w_vec <- runif(n, -2, 2)
  z_vec <- w_vec + v
  y_vec <- z_vec + v + rnorm(n, sd = 0.1)
  z <- data.frame(z = z_vec)
  w <- data.frame(w = w_vec)

  silent <- suppressWarnings(
    crsiv(
      y = y_vec,
      z = z,
      w = w,
      method = "Landweber-Fridman",
      iterate.max = 2,
      cv = "nomad",
      cv.threshold = 0,
      basis = "additive",
      degree.max = 2,
      segments.max = 2,
      nmulti = 2,
      max.bb.eval = 20,
      random.seed = 37,
      display.warnings = FALSE,
      display.nomad.progress = FALSE
    )
  )

  noisy <- suppressWarnings(
    crsiv(
      y = y_vec,
      z = z,
      w = w,
      method = "Landweber-Fridman",
      iterate.max = 2,
      cv = "nomad",
      cv.threshold = 0,
      basis = "additive",
      degree.max = 2,
      segments.max = 2,
      nmulti = 2,
      max.bb.eval = 20,
      random.seed = 37,
      display.warnings = FALSE,
      display.nomad.progress = TRUE
    )
  )

  expect_equal(as.numeric(noisy$phi), as.numeric(silent$phi), tolerance = 1e-12)
  expect_equal(as.numeric(noisy$norm.stop), as.numeric(silent$norm.stop), tolerance = 1e-12)
  expect_identical(noisy$num.iterations, silent$num.iterations)
})
