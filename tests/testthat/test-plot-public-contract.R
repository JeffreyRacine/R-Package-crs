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
    grid_control = np_grid_control(xq = 0.5),
    render_control = np_render_control(style = "band")
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
})
