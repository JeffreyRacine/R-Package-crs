test_that("regression plot payload shadows legacy mean slices", {
  set.seed(11)
  n <- 40
  d <- data.frame(
    y = rnorm(n),
    x1 = runif(n),
    x2 = runif(n),
    z = factor(sample(c("a", "b"), n, replace = TRUE))
  )
  d$y <- sin(2 * pi * d$x1) + d$x2 + as.numeric(d$z == "b") + rnorm(n, sd = 0.1)

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

  modern <- getFromNamespace(".crs_plot_payload_regression", "crs")(
    model,
    num.eval = 9,
    display.nomad.progress = FALSE
  )
  legacy <- getFromNamespace(".crs_plot_payload_regression", "crs")(
    model,
    num.eval = 9,
    legacy = TRUE,
    display.nomad.progress = FALSE
  )

  expect_s3_class(modern, "crs_plot_payload")
  expect_identical(modern$route, "crs")
  expect_identical(modern$view, "fit")
  expect_identical(modern$source, "payload")
  expect_equal(names(modern$slices), names(model$xz))
  expect_equal(names(legacy$slices), names(model$xz))

  for (nm in names(model$xz)) {
    expect_equal(nrow(modern$slices[[nm]]), nrow(legacy$slices[[nm]]))
    expect_equal(modern$slices[[nm]][[nm]], legacy$slices[[nm]][[1]],
                 tolerance = 0)
    expect_equal(modern$slices[[nm]]$fit, legacy$slices[[nm]]$mean,
                 tolerance = 1e-10)
  }
})

test_that("regression plot payload carries quantile object state", {
  set.seed(12)
  n <- 36
  d <- data.frame(y = rnorm(n), x = runif(n))
  d$y <- d$x + stats::rt(n, df = 4) * 0.05

  model <- crs(
    y ~ x,
    data = d,
    cv = "none",
    degree = 3,
    segments = 1,
    tau = 0.5,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )

  modern <- getFromNamespace(".crs_plot_payload_regression", "crs")(
    model,
    num.eval = 7,
    display.nomad.progress = FALSE
  )
  legacy <- getFromNamespace(".crs_plot_payload_regression", "crs")(
    model,
    num.eval = 7,
    legacy = TRUE,
    display.nomad.progress = FALSE
  )

  expect_equal(modern$tau, 0.5)
  expect_equal(modern$slices$x$x, legacy$slices[[1]]$x, tolerance = 0)
  expect_equal(modern$slices$x$fit, legacy$slices[[1]]$mean,
               tolerance = 1e-10)
})

test_that("regression plot payload supports modern derivatives", {
  set.seed(13)
  d <- data.frame(y = rnorm(30), x = runif(30))
  model <- crs(
    y ~ x,
    data = d,
    cv = "none",
    degree = 3,
    segments = 1,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )

  modern <- getFromNamespace(".crs_plot_payload_regression", "crs")(
    model,
    deriv = 1,
    num.eval = 7,
    display.nomad.progress = FALSE
  )
  expect_s3_class(modern, "crs_plot_payload")
  expect_identical(modern$source, "payload")
  expect_identical(modern$view, "derivative")
  expect_named(modern$slices$x, c("x", "fit"))
  expect_equal(nrow(modern$slices$x), 7)

  legacy <- getFromNamespace(".crs_plot_payload_regression", "crs")(
    model,
    deriv = 1,
    num.eval = 7,
    legacy = TRUE,
    display.nomad.progress = FALSE
  )
  expect_s3_class(legacy, "crs_plot_payload")
  expect_identical(legacy$source, "legacy")
  expect_identical(legacy$view, "derivative")
})

test_that("2D regression payload returns object-fed surface data", {
  set.seed(14)
  n <- 40
  d <- data.frame(x1 = runif(n), x2 = runif(n))
  d$y <- d$x1 - d$x2 + rnorm(n, sd = 0.05)

  model <- crs(
    y ~ x1 + x2,
    data = d,
    cv = "none",
    degree = c(3, 3),
    segments = c(1, 1),
    basis = "tensor",
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )

  payload <- getFromNamespace(".crs_plot_payload_regression", "crs")(
    model,
    num.eval = 6,
    perspective = TRUE,
    display.nomad.progress = FALSE
  )

  expect_s3_class(payload, "crs_plot_payload")
  expect_true(payload$perspective)
  expect_equal(length(payload$x), 6)
  expect_equal(length(payload$y), 6)
  expect_equal(dim(payload$z), c(6, 6))
  expect_named(payload$data[1:2], c("x1", "x2"))
  expect_equal(as.numeric(payload$z), payload$data$fit, tolerance = 0)
})

test_that("IV and IV derivative payloads expose sorted object state", {
  z <- c(0.3, 0.1, 0.2)
  iv <- list(
    num.x = 1L,
    num.z = NULL,
    xz = data.frame(z = z),
    y = c(2, 1, 3),
    phi = c(30, 10, 20),
    deriv.mat = matrix(c(3, 1, 2), ncol = 1),
    deriv.mat.lwr = matrix(c(2.5, 0.5, 1.5), ncol = 1),
    deriv.mat.upr = matrix(c(3.5, 1.5, 2.5), ncol = 1),
    xnames = "z"
  )
  class(iv) <- c("crsiv", "crs")

  fit.payload <- getFromNamespace(".crs_plot_payload_iv", "crs")(iv)
  deriv.payload <- getFromNamespace(".crs_plot_payload_iv", "crs")(
    iv,
    deriv = TRUE,
    ci = TRUE
  )

  expect_equal(fit.payload$data$z, sort(z), tolerance = 0)
  expect_equal(fit.payload$data$fit, c(10, 20, 30), tolerance = 0)
  expect_equal(deriv.payload$data$fit, c(1, 2, 3), tolerance = 0)
  expect_equal(deriv.payload$data$lwr, c(0.5, 1.5, 2.5), tolerance = 0)

  ivd <- list(
    num.x = 1L,
    num.z = NULL,
    xz = data.frame(z = z),
    phi = c(30, 10, 20),
    phi.prime = c(3, 1, 2),
    xnames = "z"
  )
  class(ivd) <- c("crsivderiv", "crs")

  payload <- getFromNamespace(".crs_plot_payload_iv_deriv", "crs")(ivd)
  phi.payload <- getFromNamespace(".crs_plot_payload_iv_deriv", "crs")(
    ivd,
    phi = TRUE
  )
  expect_equal(payload$data$fit, c(1, 2, 3), tolerance = 0)
  expect_equal(phi.payload$data$fit, c(10, 20, 30), tolerance = 0)
})

test_that("CLSD payload mirrors density, distribution, and derivative choices", {
  obj <- list(
    x = c(0.3, 0.1, 0.2),
    xer = c(0.25, 0.05, 0.15),
    density = c(3, 1, 2),
    distribution = c(0.3, 0.1, 0.2),
    density.deriv = c(30, 10, 20),
    density.er = c(4, 2, 3),
    distribution.er = c(0.4, 0.2, 0.3),
    density.deriv.er = c(40, 20, 30)
  )
  class(obj) <- "clsd"

  density <- getFromNamespace(".crs_plot_payload_clsd", "crs")(obj, er = FALSE)
  distribution <- getFromNamespace(".crs_plot_payload_clsd", "crs")(
    obj,
    distribution = TRUE
  )
  derivative <- getFromNamespace(".crs_plot_payload_clsd", "crs")(
    obj,
    derivative = TRUE
  )

  expect_equal(density$data$x, sort(obj$x), tolerance = 0)
  expect_equal(density$data$y, c(1, 2, 3), tolerance = 0)
  expect_identical(distribution$view, "distribution")
  expect_equal(distribution$data$y, c(0.2, 0.3, 0.4), tolerance = 0)
  expect_identical(derivative$view, "derivative")
  expect_equal(derivative$data$y, c(20, 30, 40), tolerance = 0)
})
