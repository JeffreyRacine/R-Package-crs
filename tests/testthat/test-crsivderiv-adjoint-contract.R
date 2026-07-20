test_that("crsivderiv Gaussian integral helper uses ordinary CDF scaling", {
  helper <- getFromNamespace(".crsiv_gaussian_integral_apply", "crs")
  z <- c(-1, 1)
  rhs <- c(2, 5)
  h <- 0.5

  observed <- helper(z.train = z, z.eval = 0, rhs = rhs, bw = h)
  expected <- drop(stats::pnorm(outer(0, z, "-") / h) %*% rhs) / length(z)

  expect_equal(observed, expected, tolerance = 1e-14)
  expect_false(isTRUE(all.equal(observed, h * expected, tolerance = 1e-14)))
})

test_that("crsivderiv Gaussian integral helper preserves matrix right-hand sides", {
  helper <- getFromNamespace(".crsiv_gaussian_integral_apply", "crs")
  z <- c(-1, 1)
  rhs <- cbind(c(2, 5), c(7, 3))
  h <- 0.5
  zeval <- c(-0.25, 0.25)

  observed <- helper(z.train = z, z.eval = zeval, rhs = rhs, bw = h)
  expected <- stats::pnorm(outer(zeval, z, "-") / h) %*% rhs / length(z)

  expect_equal(observed, expected, tolerance = 1e-14)
  expect_identical(dim(observed), c(2L, 2L))
})

test_that("crsivderiv centers both adjoint terms on the fitted residual", {
  src_path <- testthat::test_path("..", "..", "R", "crsivderiv.R")
  skip_if_not(file.exists(src_path), "source R files unavailable")
  src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")

  fitted_second_term <- gregexpr(
    "S\\.z\\*mean\\.predicted\\.model\\.E\\.mu\\.w",
    src,
    perl = TRUE
  )[[1L]]
  fitted_mean <- gregexpr(
    "mean\\.predicted\\.model\\.E\\.mu\\.w <- mean\\(predicted\\.model\\.E\\.mu\\.w\\)",
    src,
    perl = TRUE
  )[[1L]]

  expect_length(fitted_second_term[fitted_second_term > 0L], 2L)
  expect_length(fitted_mean[fitted_mean > 0L], 4L)
  expect_false(grepl("S\\.z\\*mean\\.mu", src, perl = TRUE))
  expect_false(grepl(
    "mean\\(E\\.y\\.w\\) - mean\\(predicted\\.model\\.E\\.mu\\.w\\)",
    src,
    perl = TRUE
  ))
})
