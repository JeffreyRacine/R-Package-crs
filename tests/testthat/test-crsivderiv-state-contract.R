test_that("crsivderiv state matrices contain coherent positive-index states", {
  set.seed(20260720)
  n <- 48L
  w <- runif(n, -2, 2)
  v <- rnorm(n, sd = 0.15)
  z <- 0.4 * w + v
  y <- z^2 - 0.3 * v + rnorm(n, sd = 0.05)

  fit <- suppressWarnings(
    crsivderiv(
      y = y,
      z = data.frame(z = z),
      w = data.frame(w = w),
      iterate.max = 3L,
      iterate.diff.tol = 0,
      stop.on.increase = FALSE,
      starting.values = rep(0, n),
      cv = "none",
      basis = "additive",
      display.nomad.progress = FALSE,
      display.warnings = FALSE
    )
  )

  integrate <- getFromNamespace("integrate.trapezoidal", "crs")
  expected.phi <- vapply(
    seq_len(ncol(fit$phi.prime.mat)),
    function(N) {
      value <- integrate(z, fit$phi.prime.mat[, N])
      value - mean(value) + mean(y)
    },
    numeric(n)
  )

  expect_identical(ncol(fit$phi.mat), 3L)
  expect_identical(ncol(fit$phi.prime.mat), 3L)
  expect_identical(length(fit$norm.stop), 3L)
  expect_equal(fit$phi.mat, expected.phi, tolerance = 1e-14)
  expect_equal(fit$phi, fit$phi.mat[, fit$num.iterations], tolerance = 0)
  expect_equal(fit$phi.prime,
               fit$phi.prime.mat[, fit$num.iterations], tolerance = 0)
  expect_equal(fit$starting.values.phi.prime, rep(0, n), tolerance = 0)
  expect_equal(fit$starting.values.phi, rep(mean(y), n), tolerance = 1e-14)
})

test_that("crsivderiv summary reports the selected state's criterion", {
  object <- list(
    call = quote(crsivderiv(y = y, z = z, w = w)),
    kernel = FALSE,
    tau = NULL,
    num.x = 1L,
    num.z = NULL,
    xnames = "z",
    degree = 2L,
    segments = 2L,
    include = NULL,
    lambda = NULL,
    num.x.w = NULL,
    num.z.w = NULL,
    complexity = 3L,
    knots = "uniform",
    nobs = 12L,
    num.iterations = 2L,
    norm.stop = c(0.9, 0.125, 0.8),
    nmulti = 1L,
    ptm = c(user.self = 0, sys.self = 0, elapsed = 0)
  )
  class(object) <- c("crsivderiv", "crs")

  output <- capture.output(summary(object))

  expect_true(any(grepl("Stopping rule value: 0.125", output, fixed = TRUE)))
  expect_false(any(grepl("Stopping rule value: 0.8", output, fixed = TRUE)))
})
