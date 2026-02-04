test_that("resolve_cv_maxPenalty handles cv.ls/cv.gcv/cv.aic consistently", {
  y <- 1:10
  w <- rep(1, length(y))
  n <- length(y)
  mu <- mean(y)
  mse <- mean((y - mu)^2)

  expected_ls <- 10 * (mse / (1 - 1 / n)^2)
  expected_gcv <- expected_ls
  expected_aic <- (log(mse) + (1 + 1 / n) / (1 - 3 / n)) + 10

  got_ls <- crs:::resolve_cv_maxPenalty(NULL, y, weights = w, cv.func = "cv.ls")
  got_gcv <- crs:::resolve_cv_maxPenalty(NULL, y, weights = w, cv.func = "cv.gcv")
  got_aic <- crs:::resolve_cv_maxPenalty(NULL, y, weights = w, cv.func = "cv.aic")

  expect_equal(got_ls, expected_ls, tolerance = 1e-10)
  expect_equal(got_gcv, expected_gcv, tolerance = 1e-10)
  expect_equal(got_aic, expected_aic, tolerance = 1e-10)
})

test_that("resolve_cv_maxPenalty respects weights", {
  y <- c(1, 2, 3, 4)
  w <- c(1, 2, 3, 4)
  n <- length(y)
  mu <- sum(w * y) / sum(w)
  mse <- sum(w * (y - mu)^2) / sum(w)

  expected_ls <- 10 * (mse / (1 - 1 / n)^2)
  expected_aic <- (log(mse) + (1 + 1 / n) / (1 - 3 / n)) + 10

  got_ls <- crs:::resolve_cv_maxPenalty(NULL, y, weights = w, cv.func = "cv.ls")
  got_aic <- crs:::resolve_cv_maxPenalty(NULL, y, weights = w, cv.func = "cv.aic")

  expect_equal(got_ls, expected_ls, tolerance = 1e-10)
  expect_equal(got_aic, expected_aic, tolerance = 1e-10)
})
