test_that("CRS wild bootstrap helper matches explicit hat-operator oracle", {
  set.seed(6201)
  d <- data.frame(x = runif(28), z = factor(sample(letters[1:2], 28, TRUE)))
  d$y <- sin(2 * pi * d$x) + as.numeric(d$z) / 5 + rnorm(28, sd = 0.04)
  model <- crs(
    y ~ x + z,
    data = d,
    cv = "none",
    kernel = TRUE,
    basis = "additive",
    degree = 2,
    segments = 1,
    lambda = 0.4,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )
  nd <- data.frame(
    x = seq(0.1, 0.9, length.out = 6),
    z = factor(rep(levels(d$z)[1L], 6), levels = levels(d$z))
  )

  H.train <- crshat(model, output = "matrix")
  H.eval <- crshat(model, newdata = nd, output = "matrix")
  fit.mean <- as.vector(H.train %*% model$y)
  set.seed(6202)
  helper <- getFromNamespace(".crs.bootstrap.matrix.wild", "crs")(
    object = model,
    newdata = nd,
    boot.num = 9L,
    wild = "rademacher",
    display.nomad.progress = FALSE
  )
  set.seed(6202)
  draws <- getFromNamespace(".crs_rademacher_draws", "crs")(
    n = length(model$y),
    B = 9L
  )
  ystar <- (model$y - fit.mean) * draws
  ystar <- ystar + fit.mean
  oracle <- t(H.eval %*% ystar)

  expect_equal(helper$center, as.vector(H.eval %*% model$y),
               tolerance = 1e-10)
  expect_equal(helper$boot.mat, oracle, tolerance = 1e-10)
})

test_that("CRS wild bootstrap block route reuses draws across evaluation blocks", {
  set.seed(6205)
  d <- data.frame(x = runif(26))
  d$y <- sin(2 * pi * d$x) + rnorm(26, sd = 0.04)
  model <- crs(
    y ~ x,
    data = d,
    cv = "none",
    degree = 2,
    segments = 1,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )
  nd <- data.frame(x = seq(0.05, 0.95, length.out = 11))
  H <- crshat(model, newdata = nd, output = "matrix")
  fit.mean <- as.vector(crshat(model, output = "apply"))

  old.dense <- getOption("crs.plot.wild.dense.hat.threshold.bytes")
  old.block <- getOption("crs.plot.wild.hat.block.bytes")
  on.exit(options(crs.plot.wild.dense.hat.threshold.bytes = old.dense,
                  crs.plot.wild.hat.block.bytes = old.block),
          add = TRUE)
  options(crs.plot.wild.dense.hat.threshold.bytes = 0,
          crs.plot.wild.hat.block.bytes = 8 * length(model$y) * 3)

  set.seed(6206)
  helper <- getFromNamespace(".crs.bootstrap.matrix.wild", "crs")(
    object = model,
    newdata = nd,
    boot.num = 8L,
    wild = "rademacher",
    display.nomad.progress = FALSE
  )
  set.seed(6206)
  draws <- getFromNamespace(".crs_rademacher_draws", "crs")(
    n = length(model$y),
    B = 8L
  )
  ystar <- (model$y - fit.mean) * draws
  ystar <- ystar + fit.mean
  oracle <- t(H %*% ystar)

  expect_equal(helper$boot.mat, oracle, tolerance = 1e-10)
})

test_that("plot.crs accepts explicit wild bootstrap for 1D and surface routes", {
  set.seed(6203)
  d1 <- data.frame(x = runif(30))
  d1$y <- cos(2 * pi * d1$x) + rnorm(30, sd = 0.05)
  fit1 <- crs(
    y ~ x,
    data = d1,
    cv = "none",
    degree = 2,
    segments = 1,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )
  out1 <- plot(
    fit1,
    output = "data",
    errors = "bootstrap",
    bootstrap = "wild",
    boot_control = crs_boot_control(wild = "mammen"),
    B = 7L,
    neval = 6L
  )
  expect_named(out1[[1L]], c("x", "mean", "lwr", "upr"))

  d2 <- data.frame(x1 = runif(32), x2 = runif(32))
  d2$y <- sin(2 * pi * d2$x1) + d2$x2 + rnorm(32, sd = 0.05)
  fit2 <- crs(
    y ~ x1 + x2,
    data = d2,
    cv = "none",
    basis = "tensor",
    degree = c(2, 2),
    segments = c(1, 1),
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )
  out2 <- plot(
    fit2,
    output = "data",
    perspective = TRUE,
    errors = "bootstrap",
    bootstrap = "wild",
    band = "all",
    B = 5L,
    neval = 4L
  )
  expect_named(out2[[1L]],
               c("x1", "x2", "fit", "lwr", "upr", "lwr.sim", "upr.sim",
                 "lwr.bonf", "upr.bonf"))
})

test_that("plot.crs bootstrap defaults are NP-style wild with 1999 replications", {
  fn <- getFromNamespace(".crs_plot_regression_1d_public", "crs")
  body.txt <- paste(deparse(body(fn)), collapse = "\n")
  expect_match(body.txt, 'dots\\$plot.errors.boot.method,\\s*"wild"',
               fixed = FALSE)
  expect_match(body.txt, 'dots\\$plot.errors.boot.num,\\s*1999L',
               fixed = FALSE)
})

test_that("plot.crs rejects unsupported bootstrap selectors explicitly", {
  set.seed(6204)
  d <- data.frame(x = runif(24), y = rnorm(24))
  fit <- crs(
    y ~ x,
    data = d,
    cv = "none",
    degree = 2,
    segments = 1,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )
  expect_error(plot(fit, output = "data", errors = "bootstrap",
                    bootstrap = "fixed", B = 3L),
               "bootstrap=\"wild\" or bootstrap=\"inid\"")
})
