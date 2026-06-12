test_that("plot.crs rejects stale CRS plot arguments at the public boundary", {
  set.seed(101)
  d <- data.frame(x = runif(28))
  d$y <- sin(2 * pi * d$x) + rnorm(28, sd = 0.05)
  model <- crs(
    y ~ x,
    data = d,
    cv = "none",
    degree = 2,
    segments = 1,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )

  stale <- crs_plot_contract_cases()$stale
  for (arg in stale) {
    call.args <- list(model, output = "data")
    call.args[[arg]] <- TRUE
    expect_error(do.call(plot, call.args), "unused plot argument",
                 info = sprintf("argument %s should be rejected", arg))
  }
})

test_that("plot.crs accepts canonical NP-style plot controls", {
  set.seed(102)
  d <- data.frame(x = runif(30))
  d$y <- cos(2 * pi * d$x) + rnorm(30, sd = 0.05)
  model <- crs(
    y ~ x,
    data = d,
    cv = "none",
    degree = 2,
    segments = 1,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )

  out <- plot(
    model,
    output = "data",
    neval = 6,
    data_overlay = TRUE,
    data_rug = FALSE,
    layout = "current",
    common_scale = TRUE,
    xtrim = 0.05,
    grid_control = crs_grid_control(xq = 0.5)
  )
  expect_type(out, "list")
  expect_named(out[[1L]], c("x", "mean"))

  gout <- plot(
    model,
    output = "data",
    gradients = TRUE,
    gradient_order = 1,
    neval = 6
  )
  expect_type(gout, "list")
  expect_named(gout[[1L]], c("x", "deriv"))

  bout <- plot(
    model,
    output = "data",
    errors = "bootstrap",
    band = "all",
    bootstrap = "inid",
    B = 3,
    neval = 6
  )
  expect_named(bout[[1L]],
               c("x", "mean", "lwr", "upr", "lwr.sim", "upr.sim",
                 "lwr.bonf", "upr.bonf"))
  expect_silent(crs_plot_pdf(plot(
    model,
    errors = "bootstrap",
    band = "all",
    bootstrap = "inid",
    B = 3,
    neval = 6,
    legend = FALSE
  )))
  expect_silent(crs_plot_pdf(plot(
    model,
    errors = "bootstrap",
    band = "all",
    bootstrap = "inid",
    B = 3,
    neval = 6,
    legend = list(x = "bottomleft")
  )))
})

test_that("plot route defaults render without hidden legacy bridge arguments", {
  set.seed(103)
  d1 <- data.frame(x = runif(30))
  d1$y <- sin(2 * pi * d1$x) + rnorm(30, sd = 0.05)
  fit1 <- crs(
    y ~ x,
    data = d1,
    cv = "none",
    degree = 2,
    segments = 1,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )

  d2 <- data.frame(x1 = runif(32), x2 = runif(32))
  d2$y <- sin(2 * pi * d2$x1) + d2$x2 + rnorm(32, sd = 0.05)
  fit2 <- crs(
    y ~ x1 + x2,
    data = d2,
    cv = "none",
    kernel = FALSE,
    basis = "tensor",
    degree = c(2, 2),
    segments = c(1, 1),
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )

  fitq <- crs(
    y ~ x,
    data = d1,
    cv = "none",
    degree = 2,
    segments = 1,
    tau = 0.5,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )

  expect_silent(crs_plot_pdf(plot(fit1, neval = 6)))
  expect_silent(crs_plot_pdf(plot(fit2, neval = 4, renderer = "base")))
  expect_silent(crs_plot_pdf(plot(fitq, neval = 6)))
})

test_that("plot.crs consumes NP-style data_rug and uses overlay data range", {
  set.seed(105)
  d <- data.frame(x = runif(34))
  d$y <- sin(2 * pi * d$x) + rnorm(34, sd = 0.05)
  model <- crs(
    y ~ x,
    data = d,
    cv = "none",
    degree = 2,
    segments = 1,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )
  model$y[1L] <- 100

  expect_warning(
    crs_plot_pdf(plot(model, neval = 7, data_rug = TRUE)),
    NA
  )
  overlay_range <- getFromNamespace(".crs_plot_overlay_range", "crs")
  expect_equal(overlay_range(c(-1, 1), model$y),
               c(min(-1, min(model$y)), max(model$y)))
})

test_that("plot.crs refit gradient bootstrap uses overall derivative columns", {
  set.seed(107)
  d <- data.frame(
    f = factor(rep(c("a", "b"), each = 18)),
    x = runif(36)
  )
  d$y <- 1 + 0.5 * (d$f == "b") + sin(2 * pi * d$x) + rnorm(36, sd = 0.04)
  fit <- crs(
    y ~ f + x,
    data = d,
    cv = "none",
    degree = 2,
    segments = 1,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )

  out <- suppressWarnings(plot(
    fit,
    output = "data",
    gradients = TRUE,
    gradient_order = 1,
    errors = "bootstrap",
    bootstrap = "inid",
    B = 2L,
    neval = 5L
  ))
  x.slice <- out[["x"]]
  expect_named(x.slice, c("x", "deriv", "lwr", "upr"))

  nd <- data.frame(
    f = factor(rep(getFromNamespace("uocquantile", "crs")(fit$xz$f, prob = 0.5),
                   nrow(x.slice)),
               levels = levels(fit$xz$f)),
    x = x.slice$x
  )
  fit.deriv <- fit
  fit.deriv$deriv <- 1L
  pred <- predict(fit.deriv, newdata = nd)
  x.col <- match("x", names(fit$xz))
  f.col <- match("f", names(fit$xz))

  expect_equal(x.slice$deriv, attr(pred, "deriv.mat")[, x.col],
               tolerance = 1e-10)
  expect_false(isTRUE(all.equal(x.slice$deriv,
                                attr(pred, "deriv.mat")[, f.col],
                                tolerance = 1e-7)))
})

test_that("plot.crs output=data has no graphics device side effect", {
  set.seed(108)
  d <- data.frame(x = runif(28))
  d$y <- sin(2 * pi * d$x) + rnorm(28, sd = 0.05)
  fit <- crs(
    y ~ x,
    data = d,
    cv = "none",
    degree = 2,
    segments = 1,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )

  oldwd <- getwd()
  tmp <- tempfile("crs-plot-data-")
  dir.create(tmp)
  on.exit(setwd(oldwd), add = TRUE)
  setwd(tmp)
  expect_type(plot(fit, output = "data", neval = 5), "list")
  expect_false(file.exists(file.path(tmp, "Rplots.pdf")))
})

test_that("plot.crs surface zlim follows NP-style data overlay range", {
  set.seed(106)
  d <- data.frame(x1 = runif(36), x2 = runif(36))
  d$y <- sin(2 * pi * d$x1) + d$x2 + rnorm(36, sd = 0.05)
  model <- crs(
    y ~ x1 + x2,
    data = d,
    cv = "none",
    kernel = FALSE,
    basis = "tensor",
    degree = c(2, 2),
    segments = c(1, 1),
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )
  model$y[1L] <- 100

  expect_warning(
    crs_plot_pdf(plot(model, neval = 5, perspective = TRUE,
                      data_rug = TRUE)),
    NA
  )
  payload <- getFromNamespace(".crs_plot_payload_regression", "crs")(
    object = model,
    deriv = 0L,
    ci = FALSE,
    num.eval = 5L,
    xtrim = 0,
    perspective = TRUE,
    legacy = FALSE,
    display.nomad.progress = FALSE
  )
  overlay_range <- getFromNamespace(".crs_plot_overlay_range", "crs")
  expect_equal(overlay_range(range(payload$z, finite = TRUE), model$y)[2L],
               100)
})

test_that("curve plot routes keep their documented NP-compatible controls", {
  set.seed(104)
  n <- 32
  v <- rnorm(n, sd = 0.08)
  w_vec <- runif(n, -1.5, 1.5)
  z_vec <- w_vec + v
  y_vec <- z_vec^2 + v + rnorm(n, sd = 0.10)
  fit.iv <- crsiv(
    y = y_vec,
    z = data.frame(z = z_vec),
    w = data.frame(w = w_vec),
    method = "Tikhonov",
    alpha = 0.1,
    cv = "none",
    basis = "additive",
    deriv = 1,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )
  fit.ivd <- crsivderiv(
    y = y_vec,
    z = data.frame(z = z_vec),
    w = data.frame(w = w_vec),
    iterate.max = 2,
    cv = "none",
    basis = "additive",
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )
  fit.clsd <- clsd(
    rnorm(36),
    degree = 2,
    segments = 1,
    xeval = seq(-2, 2, length.out = 9),
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )

  expect_s3_class(plot(fit.iv, output = "data"), "data.frame")
  expect_s3_class(plot(fit.iv, output = "data", deriv = TRUE), "data.frame")
  expect_s3_class(plot(fit.ivd, output = "data"), "data.frame")
  expect_s3_class(plot(fit.ivd, output = "data", phi = TRUE), "data.frame")
  expect_s3_class(plot(fit.clsd, output = "data"), "data.frame")
  expect_s3_class(plot(fit.clsd, output = "data", distribution = TRUE),
                  "data.frame")
  expect_s3_class(plot(fit.clsd, output = "data", derivative = TRUE),
                  "data.frame")

  expect_error(plot(fit.iv, output = "data", errors = "bootstrap"),
               "does not support bootstrap errors")
  expect_error(plot(fit.iv, output = "data", errors = "asymptotic"),
               "asymptotic intervals are unavailable")
  expect_error(plot(fit.iv, output = "data", data_rug = TRUE),
               "plot.crsiv does not support plot argument data_rug")
  expect_error(plot(fit.iv, output = "data", renderer = "rgl"),
               "plot.crsiv does not support plot argument renderer")
  expect_error(plot(fit.iv, output = "data", B = 5),
               "plot.crsiv does not support plot argument B")
  expect_error(plot(fit.ivd, output = "data", errors = "bootstrap"),
               "plot.crsivderiv does not support plot argument errors")
  expect_error(plot(fit.ivd, output = "data", boot_control = crs_boot_control()),
               "plot.crsivderiv does not support plot argument boot_control")
  expect_error(plot(fit.clsd, output = "data", data_overlay = TRUE),
               "plot.clsd does not support plot argument data_overlay")
  expect_error(plot(fit.clsd, output = "data", legend = TRUE),
               "plot.clsd does not support plot argument legend")
  expect_warning(
    crs_plot_pdf(plot(fit.clsd, ylab = "Density", xlab = "Support",
                      ylim = c(0, 1), type = "l")),
    NA
  )
})
