test_that("shadow 1D regression route matches legacy plot-data oracle", {
  set.seed(21)
  n <- 45
  d <- data.frame(
    x1 = runif(n),
    x2 = runif(n),
    z = factor(sample(c("a", "b"), n, replace = TRUE))
  )
  d$y <- cos(2 * pi * d$x1) + d$x2 + as.numeric(d$z == "b") +
    rnorm(n, sd = 0.1)

  model <- crs(
    y ~ x1 + x2 + z,
    data = d,
    cv = "none",
    degree = c(3, 3),
    segments = c(1, 1),
    include = 1,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )

  shadow <- getFromNamespace(".crs_plot_regression_1d_shadow", "crs")(
    model,
    output = "data",
    neval = 8
  )
  legacy <- plot(
    model,
    mean = TRUE,
    plot.behavior = "data",
    num.eval = 8,
    display.nomad.progress = FALSE
  )
  names(legacy) <- names(model$xz)

  expect_equal(names(shadow), names(legacy))
  for (nm in names(shadow)) {
    expect_equal(shadow[[nm]], legacy[[nm]], tolerance = 1e-10)
  }
})

test_that("shadow 1D regression route renders without changing public plot", {
  set.seed(22)
  d <- data.frame(x = runif(40))
  d$y <- sin(2 * pi * d$x) + rnorm(40, sd = 0.05)

  model <- crs(
    y ~ x,
    data = d,
    cv = "none",
    degree = 3,
    segments = 1,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )

  pdf.file <- tempfile(fileext = ".pdf")
  grDevices::pdf(pdf.file)
  on.exit({
    grDevices::dev.off()
    unlink(pdf.file)
  }, add = TRUE)

  expect_error(
    invisible(getFromNamespace(".crs_plot_regression_1d_shadow", "crs")(
      model,
      output = "plot",
      neval = 9,
      data_overlay = TRUE
    )),
    NA
  )
  expect_true(file.exists(pdf.file))
})

test_that("shadow 1D regression route supports interval payloads", {
  set.seed(23)
  d <- data.frame(x = runif(36))
  d$y <- d$x + rnorm(36, sd = 0.05)

  model <- crs(
    y ~ x,
    data = d,
    cv = "none",
    degree = 3,
    segments = 1,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )

  shadow <- getFromNamespace(".crs_plot_regression_1d_shadow", "crs")(
    model,
    output = "data",
    neval = 7,
    ci = TRUE
  )
  legacy <- plot(
    model,
    mean = TRUE,
    plot.behavior = "data",
    num.eval = 7,
    ci = TRUE,
    display.nomad.progress = FALSE
  )

  expect_named(shadow[[1]], c("x", "mean", "lwr", "upr"))
  expect_equal(shadow[[1]], legacy[[1]], tolerance = 1e-10)
})

test_that("shadow 1D regression route rejects unknown public dots", {
  set.seed(24)
  d <- data.frame(x = runif(20), y = rnorm(20))
  model <- crs(
    y ~ x,
    data = d,
    cv = "none",
    degree = 3,
    segments = 1,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )

  expect_error(
    getFromNamespace(".crs_plot_regression_1d_shadow", "crs")(
      model,
      output = "data",
      neva = 8
    ),
    "unused plot argument"
  )
})

test_that("public opt-in fit route matches legacy mean data oracle", {
  set.seed(25)
  n <- 42
  d <- data.frame(
    x1 = runif(n),
    x2 = runif(n),
    z = factor(sample(c("a", "b"), n, replace = TRUE))
  )
  d$y <- d$x1^2 - d$x2 + as.numeric(d$z == "b") + rnorm(n, sd = 0.05)

  model <- crs(
    y ~ x1 + x2 + z,
    data = d,
    cv = "none",
    degree = c(3, 3),
    segments = c(1, 1),
    include = 1,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )

  modern <- plot(
    model,
    plot.view = "fit",
    output = "data",
    neval = 8,
    display.nomad.progress = FALSE
  )
  legacy <- plot(
    model,
    mean = TRUE,
    plot.behavior = "data",
    num.eval = 8,
    display.nomad.progress = FALSE
  )
  names(legacy) <- names(model$xz)

  expect_equal(names(modern), names(legacy))
  for (nm in names(modern)) {
    expect_equal(modern[[nm]], legacy[[nm]], tolerance = 1e-10)
  }
})

test_that("public opt-in fit route rejects unsupported modern combinations", {
  set.seed(26)
  d <- data.frame(x = runif(24), y = rnorm(24))
  model <- crs(
    y ~ x,
    data = d,
    cv = "none",
    degree = 3,
    segments = 1,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )

  expect_error(
    plot(
      model,
      plot.view = "fit",
      output = "data",
      errors = "bootstrap",
      B = 3,
      ci = TRUE,
      display.nomad.progress = FALSE
    ),
    "asymptotic intervals only"
  )
  expect_error(
    plot(
      model,
      plot.view = "fit",
      output = "data",
      deriv = 1,
      display.nomad.progress = FALSE
    ),
    "fitted mean/quantile plots only"
  )
  expect_error(
    plot(
      model,
      plot.view = "fit",
      output = "data",
      perspective = TRUE,
      display.nomad.progress = FALSE
    ),
    "modern 2D regression plot route is not implemented yet"
  )
})
