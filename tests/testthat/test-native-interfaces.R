legacy_bs_des_c <- function(x, degree, nbreak, deriv, x.min, x.max, knots = NULL) {
  n <- length(x)
  ncol <- nbreak + degree - 1
  knots.int <- if (is.null(knots)) 0L else 1L
  knots.vec <- if (is.null(knots)) as.double(NULL) else as.double(knots)

  if (all(deriv == 0L)) {
    out <- .C("gsl_bspline",
              as.double(x),
              as.integer(n),
              as.integer(degree),
              as.integer(nbreak),
              as.double(x.min),
              as.double(x.max),
              knots.vec,
              as.integer(knots.int),
              Bx = double(n * ncol),
              PACKAGE = "crs")
  } else {
    out <- .C("gsl_bspline_deriv",
              as.double(x),
              as.integer(n),
              as.integer(degree),
              as.integer(nbreak),
              as.integer(deriv),
              as.integer(max(deriv)),
              as.double(x.min),
              as.double(x.max),
              knots.vec,
              as.integer(knots.int),
              Bx = double(n * ncol),
              PACKAGE = "crs")
  }

  matrix(out$Bx, nrow = n, ncol = ncol, byrow = TRUE)
}

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

test_that("bs.des .Call wrapper matches legacy .C results (uniform knots)", {
  set.seed(104)
  x <- runif(40)
  degree <- 3L
  nbreak <- 6L
  x.min <- min(x)
  x.max <- max(x)

  deriv0 <- rep(0L, length(x))
  new0 <- crs:::bs.des(x, degree = degree, nbreak = nbreak, deriv = deriv0, x.min = x.min, x.max = x.max, knots = NULL)
  old0 <- legacy_bs_des_c(x, degree = degree, nbreak = nbreak, deriv = deriv0, x.min = x.min, x.max = x.max, knots = NULL)
  expect_equal(new0, old0, tolerance = 1e-12)

  derivv <- sample.int(degree, length(x), replace = TRUE) - 1L
  newv <- crs:::bs.des(x, degree = degree, nbreak = nbreak, deriv = derivv, x.min = x.min, x.max = x.max, knots = NULL)
  oldv <- legacy_bs_des_c(x, degree = degree, nbreak = nbreak, deriv = derivv, x.min = x.min, x.max = x.max, knots = NULL)
  expect_equal(newv, oldv, tolerance = 1e-12)
})

test_that("bs.des .Call wrapper matches legacy .C results (quantile knots)", {
  set.seed(105)
  x <- runif(35)
  degree <- 2L
  nbreak <- 5L
  x.min <- min(x)
  x.max <- max(x)
  knots <- as.numeric(quantile(x, probs = seq(0, 1, length.out = nbreak), type = 7))

  deriv0 <- rep(0L, length(x))
  new0 <- crs:::bs.des(x, degree = degree, nbreak = nbreak, deriv = deriv0, x.min = x.min, x.max = x.max, knots = knots)
  old0 <- legacy_bs_des_c(x, degree = degree, nbreak = nbreak, deriv = deriv0, x.min = x.min, x.max = x.max, knots = knots)
  expect_equal(new0, old0, tolerance = 1e-12)
})
