manual_tensor_rows <- function(X) {
  n <- nrow(X[[1]])
  out <- lapply(seq_len(n), function(i) {
    Reduce(kronecker, lapply(X, function(M) M[i, ]))
  })
  do.call(rbind, out)
}

col_in_matrix <- function(col, M, tol = 1e-10) {
  any(colSums(abs(M - matrix(col, nrow = nrow(M), ncol = ncol(M)))) < tol)
}

test_that("tensor.prod.model.matrix matches manual row-wise Kronecker", {
  set.seed(201)
  n <- 20
  X <- list(
    cbind(1, runif(n), runif(n)),
    cbind(1, runif(n), runif(n), runif(n)),
    cbind(1, runif(n))
  )

  got <- tensor.prod.model.matrix(X)
  ref <- manual_tensor_rows(X)

  expect_true(is.matrix(got))
  expect_equal(dim(got), dim(ref))
  expect_equal(got, ref, tolerance = 1e-12)
})

test_that("glp.model.matrix is shape-stable and nested in tensor basis columns", {
  set.seed(202)
  n <- 25
  X <- list(
    cbind(1, runif(n), runif(n)),
    cbind(1, runif(n), runif(n), runif(n)),
    cbind(1, runif(n))
  )

  glp <- glp.model.matrix(X)
  tpm <- tensor.prod.model.matrix(X)

  expect_true(is.matrix(glp))
  expect_equal(nrow(glp), n)
  expect_true(ncol(glp) <= ncol(tpm))
  expect_equal(glp[, 1], rep(1, n), tolerance = 1e-12)

  # GLP terms are selected tensor-product terms for the same marginal bases.
  in_tensor <- vapply(seq_len(ncol(glp)), function(j) col_in_matrix(glp[, j], tpm), logical(1))
  expect_true(all(in_tensor))
})

test_that("single-input special case returns the original matrix", {
  set.seed(203)
  A <- cbind(1, runif(15), runif(15))

  expect_equal(tensor.prod.model.matrix(list(A)), A, tolerance = 1e-12)
  expect_equal(glp.model.matrix(list(A)), A, tolerance = 1e-12)
})
