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
