test_that("crsiv opt-in plot route returns sorted structural data", {
  set.seed(41)
  n <- 35
  v <- rnorm(n, sd = 0.1)
  w_vec <- runif(n, -2, 2)
  y_vec <- w_vec + v + rnorm(n, sd = 0.1)
  z <- data.frame(z = w_vec)
  w <- data.frame(w = w_vec)

  model <- crsiv(
    y = y_vec,
    z = z,
    w = w,
    method = "Tikhonov",
    alpha = 0.1,
    cv = "none",
    basis = "additive",
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )

  out <- plot(model,output = "data")
  ord <- order(model$xz[, 1])

  expect_s3_class(out, "data.frame")
  expect_named(out, c("z", "fit"))
  expect_equal(out$z, model$xz[, 1][ord])
  expect_equal(out$fit, model$phi[ord])
})

test_that("crsiv opt-in derivative route returns interval data", {
  set.seed(42)
  n <- 38
  v <- rnorm(n, sd = 0.1)
  w_vec <- runif(n, -2, 2)
  z_vec <- w_vec + v
  y_vec <- z_vec^2 + v + rnorm(n, sd = 0.1)
  z <- data.frame(z = z_vec)
  w <- data.frame(w = w_vec)

  model <- crsiv(
    y = y_vec,
    z = z,
    w = w,
    method = "Tikhonov",
    alpha = 0.1,
    cv = "none",
    basis = "additive",
    deriv = 1,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )

  out <- plot(model,output = "data", deriv = TRUE,
              errors = "asymptotic")

  expect_s3_class(out, "data.frame")
  expect_named(out, c("z", "fit", "lwr", "upr"))
  expect_true(all(is.finite(out$fit)))
})

test_that("crsivderiv opt-in plot route exposes derivative and phi data", {
  set.seed(43)
  n <- 42
  w_vec <- rnorm(n)
  z_vec <- 0.2 * w_vec + rnorm(n, sd = 0.1)
  y_vec <- z_vec^2 + rnorm(n, sd = 0.1)
  z <- data.frame(z = z_vec)
  w <- data.frame(w = w_vec)

  model <- crsivderiv(
    y = y_vec,
    z = z,
    w = w,
    iterate.max = 2,
    cv = "none",
    basis = "additive",
    display.nomad.progress = FALSE,
    display.warnings = FALSE
  )

  deriv.out <- plot(model,output = "data")
  phi.out <- plot(model,output = "data", phi = TRUE)
  ord <- order(model$xz[, 1])

  expect_equal(deriv.out$z, model$xz[, 1][ord])
  expect_equal(deriv.out$fit, model$phi.prime[ord])
  expect_equal(phi.out$fit, model$phi[ord])
})

test_that("clsd opt-in plot route exposes selected density-family data", {
  set.seed(44)
  x <- rnorm(80)
  model <- clsd(
    x,
    degree = 2,
    segments = 1,
    xeval = seq(-2, 2, length = 10),
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )

  dens <- plot(model,output = "data")
  dist <- plot(model,output = "data",
               distribution = TRUE)
  deriv <- plot(model,output = "data",
                derivative = TRUE)

  expect_named(dens, c("x", "fit"))
  expect_named(dist, c("x", "fit"))
  expect_named(deriv, c("x", "fit"))
  expect_true(all(diff(dens$x) >= 0))
  expect_true(all(diff(dist$x) >= 0))
  expect_true(all(diff(deriv$x) >= 0))
})

test_that("curve opt-in plot routes render to a graphics device", {
  set.seed(45)
  x <- rnorm(60)
  model <- clsd(
    x,
    degree = 2,
    segments = 1,
    xeval = seq(-2, 2, length = 8),
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
    invisible(plot(model,output = "plot")),
    NA
  )
  expect_true(file.exists(pdf.file))
})
