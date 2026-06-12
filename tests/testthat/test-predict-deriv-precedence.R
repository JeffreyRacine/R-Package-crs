test_that("predict.crs honors explicit derivative order for newdata", {
  set.seed(6201)
  d <- data.frame(x = runif(36))
  d$y <- sin(2 * pi * d$x) + rnorm(36, sd = 0.03)
  fit <- crs(
    y ~ x,
    data = d,
    cv = "none",
    degree = 3,
    segments = 1,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )
  nd <- data.frame(x = seq(0.1, 0.9, length.out = 7))

  pred.default <- predict(fit, newdata = nd)
  expect_null(attr(pred.default, "deriv.mat"))

  pred.deriv <- predict(fit, newdata = nd, deriv = 1)
  expect_equal(dim(attr(pred.deriv, "deriv.mat")), c(nrow(nd), ncol(nd)))

  fit.oracle <- fit
  fit.oracle$deriv <- 1L
  pred.oracle <- predict(fit.oracle, newdata = nd)
  expect_equal(as.numeric(pred.deriv), as.numeric(pred.oracle),
               tolerance = 1e-12)
  expect_equal(attr(pred.deriv, "deriv.mat"),
               attr(pred.oracle, "deriv.mat"),
               tolerance = 1e-12)
  expect_equal(fit$deriv, 0)
})

test_that("predict.crs explicit deriv zero suppresses stored derivatives", {
  set.seed(6202)
  d <- data.frame(x = runif(34))
  d$y <- cos(2 * pi * d$x) + rnorm(34, sd = 0.04)
  fit <- crs(
    y ~ x,
    data = d,
    cv = "none",
    degree = 3,
    segments = 1,
    deriv = 1,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )
  nd <- data.frame(x = seq(0.15, 0.85, length.out = 6))

  pred.default <- predict(fit, newdata = nd)
  expect_false(is.null(attr(pred.default, "deriv.mat")))

  pred.zero <- predict(fit, newdata = nd, deriv = 0)
  expect_null(attr(pred.zero, "deriv.mat"))
  expect_null(attr(pred.zero, "deriv.mat.lwr"))
  expect_null(attr(pred.zero, "deriv.mat.upr"))
  expect_equal(as.numeric(pred.zero), as.numeric(pred.default),
               tolerance = 1e-12)

  pred.train.zero <- predict(fit, deriv = 0)
  expect_null(attr(pred.train.zero, "deriv.mat"))
  expect_equal(as.numeric(pred.train.zero), as.numeric(fitted(fit)),
               tolerance = 1e-10)
})

test_that("predict.crs computes explicit no-newdata derivatives at training data", {
  set.seed(6203)
  d <- data.frame(x = runif(32))
  d$y <- sin(2 * pi * d$x) + rnorm(32, sd = 0.03)
  fit <- crs(
    y ~ x,
    data = d,
    cv = "none",
    degree = 3,
    segments = 1,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )

  pred.train <- predict(fit, deriv = 1)
  pred.newdata <- predict(fit, newdata = fit$xz, deriv = 1)

  expect_equal(as.numeric(pred.train), as.numeric(pred.newdata),
               tolerance = 1e-12)
  expect_equal(attr(pred.train, "deriv.mat"),
               attr(pred.newdata, "deriv.mat"),
               tolerance = 1e-12)
  expect_equal(dim(attr(pred.train, "deriv.mat")),
               c(nrow(fit$xz), ncol(fit$xz)))
})

test_that("predict.crs derivative columns preserve overall predictor order", {
  set.seed(6204)
  d <- data.frame(
    f = factor(rep(c("a", "b"), each = 20)),
    x = runif(40)
  )
  d$y <- 1 + 0.4 * (d$f == "b") + sin(2 * pi * d$x) + rnorm(40, sd = 0.03)
  fit <- crs(
    y ~ f + x,
    data = d,
    cv = "none",
    degree = 3,
    segments = 1,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )
  nd <- data.frame(
    f = factor(rep(c("a", "b"), length.out = 8), levels = levels(d$f)),
    x = seq(0.1, 0.9, length.out = 8)
  )

  pred <- predict(fit, newdata = nd, deriv = 1)
  expect_equal(colnames(nd), names(fit$xz))
  expect_equal(dim(attr(pred, "deriv.mat")), c(nrow(nd), ncol(nd)))

  fit.oracle <- fit
  fit.oracle$deriv <- 1L
  pred.oracle <- predict(fit.oracle, newdata = nd)
  expect_equal(attr(pred, "deriv.mat"),
               attr(pred.oracle, "deriv.mat"),
               tolerance = 1e-12)
  expect_false(isTRUE(all.equal(attr(pred, "deriv.mat")[, 1],
                                attr(pred, "deriv.mat")[, 2],
                                tolerance = 1e-8)))
})

test_that("predict.crs explicit derivative order works for kernel and quantile routes", {
  set.seed(6205)
  d <- data.frame(
    f = factor(sample(c("a", "b"), 38, replace = TRUE)),
    x = runif(38)
  )
  d$y <- 0.3 * (d$f == "b") + d$x^2 + rnorm(38, sd = 0.04)
  nd <- data.frame(
    f = factor(rep(c("a", "b"), length.out = 6), levels = levels(d$f)),
    x = seq(0.15, 0.85, length.out = 6)
  )

  kernel.fit <- crs(
    y ~ f + x,
    data = d,
    cv = "none",
    kernel = TRUE,
    degree = 3,
    segments = 1,
    lambda = 0.3,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )
  kernel.pred <- predict(kernel.fit, newdata = nd, deriv = 1)
  expect_equal(dim(attr(kernel.pred, "deriv.mat")), c(nrow(nd), ncol(nd)))

  q.fit <- crs(
    y ~ x,
    data = d,
    cv = "none",
    degree = 3,
    segments = 1,
    tau = 0.5,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )
  q.pred <- predict(q.fit, newdata = nd["x"], deriv = 1)
  expect_equal(dim(attr(q.pred, "deriv.mat")), c(nrow(nd), 1L))
})

test_that("predict.crs validates explicit derivative order", {
  set.seed(6206)
  d <- data.frame(x = runif(20))
  d$y <- d$x + rnorm(20, sd = 0.02)
  fit <- crs(
    y ~ x,
    data = d,
    cv = "none",
    degree = 2,
    segments = 1,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )
  nd <- data.frame(x = c(0.2, 0.4))

  expect_error(predict(fit, newdata = nd, deriv = -1),
               "non-negative integer scalar")
  expect_error(predict(fit, newdata = nd, deriv = c(0, 1)),
               "non-negative integer scalar")
  expect_error(predict(fit, newdata = nd, deriv = NA_real_),
               "non-negative integer scalar")
  expect_error(predict(fit, newdata = nd, deriv = 1.5),
               "non-negative integer scalar")
  expect_error(predict(fit, newdata = nd, deriv = "1"),
               "non-negative integer scalar")
})
