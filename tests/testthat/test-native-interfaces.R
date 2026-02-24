test_that("uniquecombs preserves row mapping via index attribute", {
  set.seed(101)
  x <- cbind(round(runif(24), 2), sample(1:4, 24, replace = TRUE))

  ux <- crs:::uniquecombs(x)
  idx <- attr(ux, "index")

  expect_true(is.matrix(ux))
  expect_length(idx, nrow(x))
  expect_true(all(idx >= 1))
  expect_true(all(idx <= nrow(ux)))
  expect_equal(unname(ux[idx, , drop = FALSE]), unname(x))
})

test_that("gsl.bs basis and predict paths return stable dimensions", {
  set.seed(102)
  x <- sort(runif(60))
  newx <- seq(min(x), max(x), length.out = 17)

  b0 <- crs::gsl.bs(x, nbreak = 7, degree = 3, deriv = 0)
  expect_s3_class(b0, "gsl.bs")
  expect_equal(nrow(b0), length(x))

  p0 <- predict(b0, newx = newx)
  expect_true(is.matrix(p0))
  expect_equal(nrow(p0), length(newx))
  expect_equal(ncol(p0), ncol(b0))

  b1 <- crs::gsl.bs(x, nbreak = 7, degree = 3, deriv = 1)
  p1 <- predict(b1, newx = newx)
  expect_true(is.matrix(p1))
  expect_equal(nrow(p1), length(newx))
  expect_equal(ncol(p1), ncol(b1))
})

test_that("hat.from.lm.fit matches qr-based hat values", {
  set.seed(103)
  n <- 80
  x1 <- runif(n)
  x2 <- runif(n)
  y <- 1 + 2 * x1 - 0.5 * x2 + rnorm(n, sd = 0.1)

  mm <- model.matrix(~ x1 + x2)
  fit <- .lm.fit(mm, y)

  qr_obj <- list(qr = fit$qr,
                 qraux = fit$qraux,
                 pivot = fit$pivot,
                 tol = fit$tol,
                 rank = fit$rank)
  class(qr_obj) <- "qr"

  h_ref <- hat(qr_obj)
  h_native <- crs:::hat.from.lm.fit(fit)

  expect_equal(h_native, h_ref, tolerance = 1e-10)
  expect_true(all(is.finite(h_native)))
})
