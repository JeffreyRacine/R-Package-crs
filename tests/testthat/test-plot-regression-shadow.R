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
    output = "data",
    neval = 8,
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
    errors = "asymptotic"
  )
  legacy <- plot(
    model,
    output = "data",
    neval = 7,
    errors = "asymptotic",
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

  forced <- FALSE
  expect_error(
    plot(
      model,
      output = "data",
      foo = {
        forced <<- TRUE
        stop("should not evaluate")
      }
    ),
    "unused plot argument: foo"
  )
  expect_false(forced)
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
    output = "data",
    neval = 8,
    display.nomad.progress = FALSE
  )
  legacy <- plot(
    model,
    output = "data",
    neval = 8,
    display.nomad.progress = FALSE
  )
  names(legacy) <- names(model$xz)

  expect_equal(names(modern), names(legacy))
  for (nm in names(modern)) {
    expect_equal(modern[[nm]], legacy[[nm]], tolerance = 1e-10)
  }
})

test_that("plot.crs defaults to fitted function data route", {
  set.seed(33)
  d <- data.frame(x = runif(34))
  d$y <- sin(2 * pi * d$x) + rnorm(34, sd = 0.05)

  model <- crs(
    y ~ x,
    data = d,
    cv = "none",
    degree = 3,
    segments = 1,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )

  default <- plot(
    model,
    output = "data",
    neval = 8,
    display.nomad.progress = FALSE
  )
  explicit <- plot(
    model,
    output = "data",
    neval = 8,
    display.nomad.progress = FALSE
  )

  expect_named(default[[1]], c("x", "mean"))
  expect_equal(default, explicit, tolerance = 1e-10)
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
    output = "data",
    neval = 7,
    errors = "asymptotic",
    display.nomad.progress = FALSE
  )
  legacy <- plot(
    model,
    output = "data",
    neval = 7,
    errors = "asymptotic",
    display.nomad.progress = FALSE
  )

  expect_equal(modern[[1]], legacy[[1]], tolerance = 1e-10)
})

test_that("public fit route supports mean bootstrap and derivative data", {
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

  boot <- suppressWarnings(plot(
    model,
    output = "data",
    errors = "bootstrap",
    B = 3,
    neval = 6,
    display.nomad.progress = FALSE
  ))
  expect_named(boot[[1]], c("x", "mean", "lwr", "upr"))
  expect_equal(nrow(boot[[1]]), 6)

  grad <- suppressWarnings(plot(
    model,
    output = "data",
    gradients = TRUE,
    neval = 6,
    display.nomad.progress = FALSE
  ))
  expect_named(grad[[1]], c("x", "deriv"))
  expect_equal(nrow(grad[[1]]), 6)

  expect_error(
    plot(
      model,
      output = "data",
      gradients = TRUE,
      errors = "bootstrap",
      B = 3,
      display.nomad.progress = FALSE
    ),
    "bootstrap intervals for derivative plots"
  )
  expect_error(
    plot(
      model,
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
    perspective = TRUE,
    output = "data",
    neval = 6,
    display.nomad.progress = FALSE
  )
  legacy <- plot(
    model,
    perspective = TRUE, renderer = "rgl",
    output = "data",
    neval = 6,
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
      perspective = TRUE,
      renderer = "base",
      output = "data",
      display.nomad.progress = FALSE
    ),
    NA
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
    output = "data",
    neval = 7,
    display.nomad.progress = FALSE
  ))
  legacy <- suppressWarnings(plot(
    model,
    output = "data",
    neval = 7,
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
    perspective = TRUE,
    output = "data",
    neval = 6,
    display.nomad.progress = FALSE
  ))
  legacy <- suppressWarnings(plot(
    model,
    perspective = TRUE, renderer = "rgl",
    output = "data",
    neval = 6,
    display.nomad.progress = FALSE
  ))

  expect_equal(modern, legacy, tolerance = 1e-10)
})
