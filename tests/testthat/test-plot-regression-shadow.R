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

test_that("public opt-in fit route matches legacy asymptotic interval data", {
  set.seed(32)
  d <- data.frame(x = runif(38))
  d$y <- 1 + d$x + rnorm(38, sd = 0.05)

  model <- crs(
    y ~ x,
    data = d,
    cv = "none",
    degree = 3,
    segments = 1,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )

  modern <- plot(
    model,
    plot.view = "fit",
    output = "data",
    neval = 7,
    ci = TRUE,
    display.nomad.progress = FALSE
  )
  legacy <- plot(
    model,
    mean = TRUE,
    plot.behavior = "data",
    num.eval = 7,
    ci = TRUE,
    display.nomad.progress = FALSE
  )

  expect_equal(modern[[1]], legacy[[1]], tolerance = 1e-10)
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
    "two continuous predictors"
  )
})

test_that("public opt-in surface route matches legacy surface data oracle", {
  set.seed(27)
  n <- 44
  d <- data.frame(x1 = runif(n), x2 = runif(n))
  d$y <- sin(2 * pi * d$x1) + d$x2 + rnorm(n, sd = 0.05)

  model <- crs(
    y ~ x1 + x2,
    data = d,
    cv = "none",
    basis = "tensor",
    degree = c(3, 3),
    segments = c(1, 1),
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )

  modern <- plot(
    model,
    plot.view = "fit",
    perspective = TRUE,
    output = "data",
    neval = 6,
    display.nomad.progress = FALSE
  )
  legacy <- plot(
    model,
    mean = TRUE,
    persp.rgl = TRUE,
    plot.behavior = "data",
    num.eval = 6,
    display.nomad.progress = FALSE
  )

  expect_equal(modern, legacy, tolerance = 1e-10)
})

test_that("public opt-in surface route renders with base persp", {
  set.seed(28)
  d <- data.frame(x1 = runif(36), x2 = runif(36))
  d$y <- d$x1 - d$x2 + rnorm(36, sd = 0.05)

  model <- crs(
    y ~ x1 + x2,
    data = d,
    cv = "none",
    basis = "tensor",
    degree = c(3, 3),
    segments = c(1, 1),
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
    invisible(plot(
      model,
      plot.view = "fit",
      perspective = TRUE,
      output = "plot",
      neval = 6,
      data_overlay = TRUE,
      display.nomad.progress = FALSE
    )),
    NA
  )
  expect_true(file.exists(pdf.file))
})

test_that("public opt-in surface route rejects contradictory renderer controls", {
  set.seed(29)
  d <- data.frame(x1 = runif(30), x2 = runif(30), y = rnorm(30))
  model <- crs(
    y ~ x1 + x2,
    data = d,
    cv = "none",
    basis = "tensor",
    degree = c(3, 3),
    segments = c(1, 1),
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )

  expect_error(
    plot(
      model,
      plot.view = "fit",
      persp.rgl = TRUE,
      renderer = "base",
      output = "data",
      display.nomad.progress = FALSE
    ),
    "cannot supply persp.rgl=TRUE"
  )
})

test_that("public opt-in fit route mirrors legacy quantile plot data", {
  set.seed(30)
  n <- 48
  d <- data.frame(x = runif(n), z = factor(sample(c("a", "b"), n, TRUE)))
  d$y <- 0.5 + d$x + as.numeric(d$z == "b") + rt(n, df = 4) * 0.1

  model <- suppressWarnings(crs(
    y ~ x + z,
    data = d,
    tau = 0.5,
    cv = "none",
    degree = 3,
    segments = 1,
    include = 1,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  ))

  modern <- suppressWarnings(plot(
    model,
    plot.view = "fit",
    output = "data",
    neval = 7,
    display.nomad.progress = FALSE
  ))
  legacy <- suppressWarnings(plot(
    model,
    mean = TRUE,
    plot.behavior = "data",
    num.eval = 7,
    display.nomad.progress = FALSE
  ))
  names(legacy) <- names(model$xz)

  expect_equal(modern, legacy, tolerance = 1e-10)
})

test_that("public opt-in surface route mirrors legacy quantile surface data", {
  set.seed(31)
  n <- 46
  d <- data.frame(x1 = runif(n), x2 = runif(n))
  d$y <- d$x1 - 0.5 * d$x2 + rt(n, df = 5) * 0.1

  model <- suppressWarnings(crs(
    y ~ x1 + x2,
    data = d,
    tau = 0.5,
    cv = "none",
    basis = "additive",
    degree = c(3, 3),
    segments = c(1, 1),
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  ))

  modern <- suppressWarnings(plot(
    model,
    plot.view = "fit",
    perspective = TRUE,
    output = "data",
    neval = 6,
    display.nomad.progress = FALSE
  ))
  legacy <- suppressWarnings(plot(
    model,
    mean = TRUE,
    persp.rgl = TRUE,
    plot.behavior = "data",
    num.eval = 6,
    display.nomad.progress = FALSE
  ))

  expect_equal(modern, legacy, tolerance = 1e-10)
})
