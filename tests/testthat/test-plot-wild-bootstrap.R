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

test_that("CRS wild bootstrap helper matches derivative and effect hat oracles", {
  set.seed(6203)
  n <- 30
  d <- data.frame(
    x = runif(n),
    z = runif(n),
    f = factor(sample(c("a", "b"), n, replace = TRUE))
  )
  d$y <- sin(2 * pi * d$x) + 0.5 * d$z + 0.3 * (d$f == "b") +
    rnorm(n, sd = 0.04)
  model <- crs(
    y ~ x + z + f,
    data = d,
    cv = "none",
    kernel = FALSE,
    basis = "additive",
    degree = c(3, 2),
    segments = c(1, 1),
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )
  nd <- data.frame(
    x = seq(0.15, 0.85, length.out = 7),
    z = mean(d$z),
    f = factor(levels(d$f)[1L], levels = levels(d$f))
  )
  H.deriv <- crshat(model, newdata = nd, output = "matrix",
                   deriv = 1, deriv.index = 1)
  fit.mean <- as.vector(crshat(model, output = "apply"))

  set.seed(6204)
  helper.deriv <- getFromNamespace(".crs.bootstrap.matrix.wild", "crs")(
    object = model,
    newdata = nd,
    deriv = 1,
    deriv.index = 1,
    boot.num = 6L,
    wild = "rademacher",
    display.nomad.progress = FALSE
  )
  set.seed(6204)
  draws <- getFromNamespace(".crs_rademacher_draws", "crs")(
    n = length(model$y),
    B = 6L
  )
  ystar <- (model$y - fit.mean) * draws
  ystar <- ystar + fit.mean
  oracle.deriv <- t(H.deriv %*% ystar)

  expect_equal(helper.deriv$center, as.vector(H.deriv %*% model$y),
               tolerance = 1e-10)
  expect_equal(helper.deriv$boot.mat, oracle.deriv, tolerance = 1e-10)

  nd.focal <- nd
  nd.base <- nd
  nd.focal$f <- factor(levels(d$f)[2L], levels = levels(d$f))
  nd.base$f <- factor(levels(d$f)[1L], levels = levels(d$f))
  H.effect <- crshat(model, newdata = nd.focal, output = "matrix") -
    crshat(model, newdata = nd.base, output = "matrix")

  set.seed(6207)
  helper.effect <- getFromNamespace(".crs.bootstrap.matrix.wild", "crs")(
    object = model,
    newdata = nd.focal,
    newdata.base = nd.base,
    boot.num = 6L,
    wild = "rademacher",
    display.nomad.progress = FALSE
  )
  set.seed(6207)
  draws.effect <- getFromNamespace(".crs_rademacher_draws", "crs")(
    n = length(model$y),
    B = 6L
  )
  ystar.effect <- (model$y - fit.mean) * draws.effect
  ystar.effect <- ystar.effect + fit.mean
  oracle.effect <- t(H.effect %*% ystar.effect)

  expect_equal(helper.effect$center, as.vector(H.effect %*% model$y),
               tolerance = 1e-10)
  expect_equal(helper.effect$boot.mat, oracle.effect, tolerance = 1e-10,
               ignore_attr = TRUE)
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

test_that("CRS bootstrap interval summary matches separate band construction", {
  set.seed(6210)
  boot.t <- matrix(rnorm(17 * 9), nrow = 17, ncol = 9)
  t0 <- colMeans(boot.t)
  summary <- getFromNamespace(".crs_plot_bootstrap_interval_summary", "crs")(
    boot.t = boot.t,
    t0 = t0,
    alpha = 0.1,
    band.type = "all",
    display.nomad.progress = FALSE
  )
  separate.pointwise <- getFromNamespace(".crs.bootstrap.bounds", "crs")(
    boot.t, 0.1, "pointwise", t0
  )
  separate.sim <- getFromNamespace(".crs.bootstrap.bounds", "crs")(
    boot.t, 0.1, "simultaneous", t0
  )
  separate.bonf <- getFromNamespace(".crs.bootstrap.bounds", "crs")(
    boot.t, 0.1, "bonferroni", t0
  )

  expect_named(summary, c("bounds", "all.bounds", "err", "all.err"))
  expect_equal(summary$bounds, summary$all.bounds$pointwise)
  expect_equal(summary$all.bounds$pointwise, separate.pointwise)
  expect_equal(summary$all.bounds$simultaneous, separate.sim)
  expect_equal(summary$all.bounds$bonferroni, separate.bonf)
  expect_equal(summary$err, cbind(t0 - summary$bounds[, 1L],
                                  summary$bounds[, 2L] - t0))
})

test_that("CRS wild bootstrap block progress reports B and interval construction", {
  test_time_values <- function(values) {
    i <- 0L
    force(values)
    function() {
      i <<- min(i + 1L, length(values))
      values[[i]]
    }
  }
  set.seed(6211)
  old_opts <- options(
    crs.messages = TRUE,
    crs.plot.wild.dense.hat.threshold.bytes = 0,
    crs.plot.wild.hat.block.bytes = 8 * 30 * 3,
    crs.plot.progress.start.grace.sec = 0,
    crs.plot.progress.interval.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  d <- data.frame(x = runif(30))
  d$y <- cos(2 * pi * d$x) + rnorm(30, sd = 0.04)
  fit <- crs(
    y ~ x,
    data = d,
    cv = "none",
    degree = 2,
    segments = 1,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )

  traced <- capture_crs_progress_shadow_trace(
    suppressWarnings(
      plot(
        fit,
        output = "data",
        errors = "bootstrap",
        bootstrap = "wild",
        band = "all",
        B = 8L,
        neval = 11L
      )
    ),
    force_renderer = "legacy",
    now = test_time_values(seq(0, 20, by = 0.05)),
    interactive = TRUE
  )

  lines <- vapply(traced$trace, `[[`, character(1L), "line")
  expect_true(any(grepl("Plot bootstrap", lines, fixed = TRUE)))
  expect_true(any(grepl("8/8", lines, fixed = TRUE)))
  expect_false(any(grepl("/4", lines, fixed = TRUE)))
  expect_true(any(grepl("Constructing bootstrap all bands", lines, fixed = TRUE)))
  expect_true(any(grepl("33/33", lines, fixed = TRUE)))
  expect_type(traced$value, "list")
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

test_that("plot.crs bootstrap defaults are object-aware with 1999 replications", {
  fn <- getFromNamespace(".crs_plot_regression_1d_public", "crs")
  body.txt <- paste(deparse(body(fn)), collapse = "\n")
  expect_match(body.txt, 'dots\\$plot.errors.boot.num,\\s*1999L',
               fixed = FALSE)

  default_method <- getFromNamespace(".crs_plot_default_bootstrap_method",
                                    "crs")
  mean.object <- list(tau = NULL)
  quantile.object <- list(tau = 0.5)
  expect_identical(default_method(mean.object), "wild")
  expect_identical(default_method(quantile.object), "inid")
})

test_that("CRS block count drawer produces valid fixed and geometric counts", {
  set.seed(6215)
  fixed.drawer <- getFromNamespace(".crs_block_counts_drawer", "crs")(
    n = 12L,
    B = 5L,
    blocklen = 3L,
    sim = "fixed"
  )
  geom.drawer <- getFromNamespace(".crs_block_counts_drawer", "crs")(
    n = 12L,
    B = 5L,
    blocklen = 3L,
    sim = "geom"
  )
  fixed <- fixed.drawer(1L, 5L)
  geom <- geom.drawer(1L, 5L)

  expect_equal(dim(fixed), c(12L, 5L))
  expect_equal(dim(geom), c(12L, 5L))
  expect_true(all(fixed >= 0))
  expect_true(all(geom >= 0))
  expect_equal(colSums(fixed), rep(12, 5))
  expect_equal(colSums(geom), rep(12, 5))
})

test_that("plot.crs accepts fixed and geom bootstraps for mean routes", {
  set.seed(6204)
  d <- data.frame(x = runif(28), y = rnorm(28))
  fit1 <- crs(
    y ~ x,
    data = d,
    cv = "none",
    degree = 2,
    segments = 1,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )
  fixed <- plot(
    fit1,
    output = "data",
    errors = "bootstrap",
    bootstrap = "fixed",
    boot_control = crs_boot_control(blocklen = 3),
    B = 4L,
    neval = 5L
  )
  geom <- plot(
    fit1,
    output = "data",
    errors = "bootstrap",
    bootstrap = "geom",
    boot_control = crs_boot_control(blocklen = 3),
    B = 4L,
    neval = 5L
  )
  expect_named(fixed[[1L]], c("x", "mean", "lwr", "upr"))
  expect_named(geom[[1L]], c("x", "mean", "lwr", "upr"))

  d2 <- data.frame(x1 = runif(30), x2 = runif(30))
  d2$y <- d2$x1 - d2$x2 + rnorm(30, sd = 0.04)
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
  surface <- plot(
    fit2,
    output = "data",
    perspective = TRUE,
    errors = "bootstrap",
    bootstrap = "fixed",
    boot_control = crs_boot_control(blocklen = 2),
    B = 3L,
    neval = 4L
  )
  expect_named(surface[[1L]], c("x1", "x2", "fit", "lwr", "upr"))
})

test_that("plot.crs consumes fixed bootstrap controls before base graphics", {
  set.seed(6205)
  d <- data.frame(x = runif(30))
  d$y <- sin(2 * pi * d$x) + rnorm(30, sd = 0.05)
  fit <- crs(
    y ~ x,
    data = d,
    cv = "none",
    degree = 2,
    segments = 1,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )

  expect_warning(
    crs_plot_pdf(plot(
      fit,
      errors = "bootstrap",
      bootstrap = "fixed",
      boot_control = crs_boot_control(blocklen = 5),
      B = 4L,
      neval = 5L
    )),
    NA
  )
})

test_that("plot.crs quantile bootstrap defaults to NPQREG-style refit selectors", {
  set.seed(6217)
  d <- data.frame(x = runif(30))
  d$y <- sin(2 * pi * d$x) + rnorm(30, sd = 0.08)
  fit <- crs(
    y ~ x,
    data = d,
    cv = "none",
    degree = 2,
    segments = 1,
    tau = 0.5,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )
  default <- suppressWarnings(plot(
    fit,
    output = "data",
    errors = "bootstrap",
    B = 3L,
    neval = 5L
  ))
  fixed <- suppressWarnings(plot(
    fit,
    output = "data",
    errors = "bootstrap",
    bootstrap = "fixed",
    boot_control = crs_boot_control(blocklen = 2),
    B = 3L,
    neval = 5L
  ))
  geom <- suppressWarnings(plot(
    fit,
    output = "data",
    errors = "bootstrap",
    bootstrap = "geom",
    boot_control = crs_boot_control(blocklen = 2),
    B = 3L,
    neval = 5L
  ))

  expect_named(default[[1L]], c("x", "mean", "lwr", "upr"))
  expect_named(fixed[[1L]], c("x", "mean", "lwr", "upr"))
  expect_named(geom[[1L]], c("x", "mean", "lwr", "upr"))
  expect_error(
    plot(
      fit,
      output = "data",
      errors = "bootstrap",
      bootstrap = "wild",
      B = 3L,
      neval = 5L
    ),
    "mean CRS objects only"
  )
})

test_that("plot.crs fixed bootstrap computes an NP-style default block length", {
  set.seed(6216)
  d <- data.frame(x = runif(26), y = rnorm(26))
  fit <- crs(
    y ~ x,
    data = d,
    cv = "none",
    degree = 2,
    segments = 1,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )
  out <- plot(
    fit,
    output = "data",
    errors = "bootstrap",
    bootstrap = "fixed",
    B = 3L,
    neval = 4L
  )
  expect_named(out[[1L]], c("x", "mean", "lwr", "upr"))
})
