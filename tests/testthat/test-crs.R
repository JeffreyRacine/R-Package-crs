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

test_that("crs default method routes nmulti through model selection", {
  set.seed(21)
  n <- 80
  x1 <- runif(n)
  x2 <- runif(n)
  z <- factor(rbinom(n, 1, .5))
  y <- cos(2 * pi * x1) + sin(2 * pi * x2) + as.numeric(z) - 1 + rnorm(n, sd = .25)
  xz <- data.frame(x1 = x1, x2 = x2, z = z)

  model <- crs(
    xz = xz,
    y = y,
    kernel = TRUE,
    nmulti = 2,
    degree.max = 2,
    segments.max = 2,
    max.bb.eval = 40,
    random.seed = 7,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )

  expect_s3_class(model, "crs")
  expect_identical(model$nmulti, 2)
  expect_equal(length(predict(model)), n)
})

test_that("crs cv.threshold and exhaustive search semantics are explicit", {
  set.seed(17)
  n <- 80
  x <- runif(n)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.2)

  default_auto <- crs(
    y ~ x,
    cv = "nomad",
    degree.max = 3,
    segments.max = 3,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )

  threshold_zero <- crs(
    y ~ x,
    cv = "nomad",
    cv.threshold = 0,
    degree.max = 3,
    segments.max = 3,
    max.bb.eval = 20,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )

  exhaustive <- crs(
    y ~ x,
    cv = "exhaustive",
    nmulti = 7,
    degree.max = 3,
    segments.max = 3,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )

  expect_identical(default_auto$cv, "exhaustive")
  expect_null(default_auto$nomad.restart.contract)

  expect_identical(threshold_zero$cv, "nomad")
  expect_false(is.null(threshold_zero$nomad.restart.contract))

  expect_identical(exhaustive$cv, "exhaustive")
  expect_null(exhaustive$nomad.restart.contract)

  exhaustive_summary <- capture.output(summary(exhaustive))
  nomad_summary <- capture.output(summary(threshold_zero))

  expect_true(any(grepl("Search method: exhaustive", exhaustive_summary, fixed = TRUE)))
  expect_false(any(grepl("Number of multistarts", exhaustive_summary, fixed = TRUE)))
  expect_true(any(grepl("Search method: nomad", nomad_summary, fixed = TRUE)))
  expect_true(any(grepl("Number of multistarts", nomad_summary, fixed = TRUE)))
})

test_that("crs NOMAD progress uses rich managed lines for kernel and factor routes", {
  test_time_values <- function(values) {
    i <- 0L
    force(values)
    function() {
      i <<- min(i + 1L, length(values))
      values[[i]]
    }
  }

  set.seed(11)
  n <- 40L
  x <- runif(n)
  z <- factor(sample(c("a", "b"), n, replace = TRUE))
  y <- sin(2 * pi * x) + as.numeric(z) - 1 + rnorm(n, sd = 0.1)
  dat <- data.frame(y = y, x = x, z = z)

  capture_lines <- function(kernel) {
    traced <- capture_crs_progress_shadow_trace(
      suppressWarnings(
        crs(
          y ~ x + z,
          data = dat,
          kernel = kernel,
          basis = "additive",
          cv = "nomad",
          degree.max = 2,
          segments.max = 2,
          nmulti = 2,
          max.bb.eval = 20,
          display.warnings = FALSE,
          display.nomad.progress = TRUE
        )
      ),
      force_renderer = "legacy",
      now = test_time_values(seq(0, 20, by = 0.25)),
      interactive = TRUE
    )

    unique(vapply(traced$trace, `[[`, character(1L), "line"))
  }

  kernel.lines <- capture_lines(TRUE)
  expect_true(any(grepl("^\\[crs\\] Selecting spline model\\.\\.\\.", kernel.lines)))
  expect_true(any(grepl("multistart 1/2", kernel.lines, fixed = TRUE)))
  expect_true(any(grepl("eval ", kernel.lines, fixed = TRUE)))
  expect_true(any(grepl("deg (", kernel.lines, fixed = TRUE)))
  expect_true(any(grepl("seg (", kernel.lines, fixed = TRUE)))
  expect_true(any(grepl("lambda (", kernel.lines, fixed = TRUE)))
  expect_true(any(grepl("best deg (", kernel.lines, fixed = TRUE)))
  expect_true(any(grepl("fv=", kernel.lines, fixed = TRUE)))
  expect_false(any(grepl("Calling NOMAD", kernel.lines, fixed = TRUE)))
  expect_false(any(grepl("NOMAD search", kernel.lines, fixed = TRUE)))

  factor.lines <- capture_lines(FALSE)
  expect_true(any(grepl("^\\[crs\\] Selecting spline model\\.\\.\\.", factor.lines)))
  expect_true(any(grepl("multistart 1/2", factor.lines, fixed = TRUE)))
  expect_true(any(grepl("eval ", factor.lines, fixed = TRUE)))
  expect_true(any(grepl("deg (", factor.lines, fixed = TRUE)))
  expect_true(any(grepl("seg (", factor.lines, fixed = TRUE)))
  expect_true(any(grepl("include (", factor.lines, fixed = TRUE)))
  expect_true(any(grepl("best deg (", factor.lines, fixed = TRUE)))
  expect_true(any(grepl("fv=", factor.lines, fixed = TRUE)))
  expect_false(any(grepl("Calling NOMAD", factor.lines, fixed = TRUE)))
  expect_false(any(grepl("NOMAD search", factor.lines, fixed = TRUE)))
})

test_that("crs NOMAD progress does not change fixed-seed results", {
  set.seed(13)
  n <- 40L
  x <- runif(n)
  z <- factor(sample(c("a", "b"), n, replace = TRUE))
  y <- sin(2 * pi * x) + as.numeric(z) - 1 + rnorm(n, sd = 0.1)
  dat <- data.frame(y = y, x = x, z = z)

  kernel.silent <- suppressWarnings(
    crs(
      y ~ x + z,
      data = dat,
      kernel = TRUE,
      basis = "additive",
      cv = "nomad",
      degree.max = 2,
      segments.max = 2,
      nmulti = 2,
      max.bb.eval = 20,
      random.seed = 29,
      display.warnings = FALSE,
      display.nomad.progress = FALSE
    )
  )

  kernel.noisy <- suppressWarnings(
    crs(
      y ~ x + z,
      data = dat,
      kernel = TRUE,
      basis = "additive",
      cv = "nomad",
      degree.max = 2,
      segments.max = 2,
      nmulti = 2,
      max.bb.eval = 20,
      random.seed = 29,
      display.warnings = FALSE,
      display.nomad.progress = TRUE
    )
  )

  expect_equal(as.numeric(kernel.noisy$degree), as.numeric(kernel.silent$degree))
  expect_equal(as.numeric(kernel.noisy$segments), as.numeric(kernel.silent$segments))
  expect_equal(as.numeric(kernel.noisy$lambda), as.numeric(kernel.silent$lambda), tolerance = 1e-12)
  expect_equal(as.numeric(kernel.noisy$cv.min), as.numeric(kernel.silent$cv.min), tolerance = 1e-12)

  factor.silent <- suppressWarnings(
    crs(
      y ~ x + z,
      data = dat,
      kernel = FALSE,
      basis = "additive",
      cv = "nomad",
      degree.max = 2,
      segments.max = 2,
      nmulti = 2,
      max.bb.eval = 20,
      random.seed = 31,
      display.warnings = FALSE,
      display.nomad.progress = FALSE
    )
  )

  factor.noisy <- suppressWarnings(
    crs(
      y ~ x + z,
      data = dat,
      kernel = FALSE,
      basis = "additive",
      cv = "nomad",
      degree.max = 2,
      segments.max = 2,
      nmulti = 2,
      max.bb.eval = 20,
      random.seed = 31,
      display.warnings = FALSE,
      display.nomad.progress = TRUE
    )
  )

  expect_equal(as.numeric(factor.noisy$degree), as.numeric(factor.silent$degree))
  expect_equal(as.numeric(factor.noisy$segments), as.numeric(factor.silent$segments))
  expect_equal(as.numeric(factor.noisy$include), as.numeric(factor.silent$include))
  expect_equal(as.numeric(factor.noisy$cv.min), as.numeric(factor.silent$cv.min), tolerance = 1e-12)
})

test_that("plot.crs bootstrap progress stays visible", {
  test_time_values <- function(values) {
    i <- 0L
    force(values)
    function() {
      i <<- min(i + 1L, length(values))
      values[[i]]
    }
  }
  set.seed(7)
  old_opts <- options(
    crs.progress.start.grace.known.sec = 0,
    crs.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)
  n <- 40
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(2 * pi * x1) + x2 + rnorm(n, sd = 0.1)

  model <- crs(
    y ~ x1 + x2,
    degree = c(3, 3),
    segments = c(1, 1),
    cv = "none",
    basis = "additive",
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )

  traced <- capture_crs_progress_shadow_trace(
    suppressWarnings(
      plot(
        model,
        errors = "bootstrap",
        B = 3,
        output = "data",
        perspective = FALSE,
        neval = 8,
        display.warnings = FALSE,
        display.nomad.progress = TRUE
      )
    ),
    force_renderer = "legacy",
    now = test_time_values(seq(0, 20, by = 0.25)),
    interactive = TRUE
  )

  lines <- vapply(traced$trace, `[[`, character(1L), "line")
  expect_true(any(grepl("Plot bootstrap", lines, fixed = TRUE)))
  expect_true(any(grepl("3/3", lines, fixed = TRUE)))
  expect_type(traced$value, "list")
})
