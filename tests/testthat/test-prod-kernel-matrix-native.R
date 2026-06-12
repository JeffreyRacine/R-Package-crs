prod_kernel_matrix_ref <- function(Z, z, lambda, is.ordered.z) {
  if (!is.matrix(Z)) Z <- as.matrix(Z)
  n <- NROW(Z)
  prodker <- rep.int(1, n)
  z_vec <- as.vector(z)
  lambda_vec <- as.vector(lambda)

  for (j in seq_len(NCOL(Z))) {
    if (!is.ordered.z[j]) {
      tmp <- rep.int(lambda_vec[j], n)
      eq <- Z[, j] == z_vec[j]
      tmp[eq] <- 1
    } else {
      tmp <- lambda_vec[j]^abs(Z[, j] - z_vec[j])
      eq <- Z[, j] == z_vec[j]
      tmp[eq] <- 1
    }
    prodker <- prodker * tmp
  }

  prodker
}

test_that("native prod.kernel.matrix matches R reference", {
  set.seed(20260611)
  grids <- list(
    list(
      Z = cbind(sample(1:4, 80, TRUE), sample(1:5, 80, TRUE)),
      z = c(2, 4),
      lambda = c(0.25, 0.8),
      is.ordered.z = c(FALSE, TRUE)
    ),
    list(
      Z = matrix(sample(c(1, 2, 3, NA), 90, TRUE), ncol = 3),
      z = c(NA, 2, 1),
      lambda = c(0.1, 0.5, 1),
      is.ordered.z = c(FALSE, FALSE, TRUE)
    ),
    list(
      Z = matrix(as.numeric(sample(0:7, 120, TRUE)), ncol = 4),
      z = matrix(c(0, 3, 7, 2), ncol = 1),
      lambda = matrix(c(0, 0.2, 0.7, 1), ncol = 1),
      is.ordered.z = c(FALSE, TRUE, TRUE, FALSE)
    )
  )

  for (case in grids) {
    expect_equal(
      crs:::prod.kernel.matrix(
        Z = case$Z,
        z = case$z,
        lambda = case$lambda,
        is.ordered.z = case$is.ordered.z
      ),
      prod_kernel_matrix_ref(
        Z = case$Z,
        z = case$z,
        lambda = case$lambda,
        is.ordered.z = case$is.ordered.z
      ),
      tolerance = 0
    )
  }
})

test_that("native prod.kernel.matrix preserves validation errors", {
  Z <- matrix(1:6, ncol = 2)
  expect_error(crs:::prod.kernel.matrix(Z, z = 1, lambda = c(0.5, 0.5), is.ordered.z = c(FALSE, TRUE)),
               "incompatible dimensions")
  expect_error(crs:::prod.kernel.matrix(Z, z = c(1, 2), lambda = c(0.5, 0.5), is.ordered.z = FALSE),
               "is.ordered.z and Z incompatible")
})

test_that("native prod.kernel.matrix pins NA behavior", {
  Z <- cbind(
    unordered = c(NA_real_, 1, 2),
    ordered = c(1, NA_real_, 2)
  )

  out <- crs:::prod.kernel.matrix(
    Z = Z,
    z = c(1, 1),
    lambda = c(0.25, 0.5),
    is.ordered.z = c(FALSE, TRUE)
  )

  expect_equal(out[1], 0.25, tolerance = 0)
  expect_true(is.na(out[2]))
  expect_equal(out[3], 0.125, tolerance = 0)
})
