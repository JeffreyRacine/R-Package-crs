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
        mean = TRUE,
        ci = TRUE,
        plot.errors.method = "bootstrap",
        plot.errors.boot.num = 3,
        plot.behavior = "data",
        num.eval = 8,
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
