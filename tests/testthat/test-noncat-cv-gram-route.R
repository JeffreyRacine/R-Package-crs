make_noncat_cv_gram_data <- function(n = 140L, seed = 20260611L) {
  set.seed(seed)
  x <- cbind(x1 = runif(n), x2 = runif(n))
  y <- sin(2 * pi * x[, 1L]) + 0.5 * cos(2 * pi * x[, 2L]) +
    rnorm(n, sd = 0.2)
  list(x = x, y = y, weights = seq(0.6, 1.4, length.out = n))
}

noncat_cv_gram_value <- function(data,
                                 basis,
                                 cv_func,
                                 weighted,
                                 use_gram) {
  K <- matrix(c(3, 2, 2, 2), ncol = 2L, byrow = TRUE)
  crs:::cv.kernel.spline(
    x = data$x,
    y = data$y,
    z = NULL,
    K = K,
    lambda = NULL,
    z.unique = NULL,
    ind = NULL,
    ind.vals = NULL,
    ind.list = NULL,
    nrow.z.unique = 0L,
    is.ordered.z = FALSE,
    knots = "quantiles",
    basis = basis,
    cv.func = cv_func,
    weights = if(weighted) data$weights else NULL,
    tau = NULL,
    singular.ok = FALSE,
    display.warnings = FALSE,
    use.ridge = FALSE,
    smooth.penalty = TRUE,
    penalty.scale = 1000,
    use.gram.cv = use_gram,
    gram.rcond.min = 1e-10,
    record.gram.stats = TRUE
  )
}

test_that("non-categorical Gram CV route preserves QR objectives", {
  cases <- expand.grid(
    basis = c("additive", "tensor", "glp"),
    cv_func = c("cv.ls", "cv.gcv", "cv.aic"),
    weighted = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )

  for(i in seq_len(nrow(cases))) {
    row <- cases[i, ]
    data <- make_noncat_cv_gram_data(seed = 91000L + i)
    direct <- noncat_cv_gram_value(data, row$basis, row$cv_func,
                                   row$weighted, use_gram = FALSE)
    gram <- noncat_cv_gram_value(data, row$basis, row$cv_func,
                                 row$weighted, use_gram = TRUE)

    expect_equal(as.numeric(gram), as.numeric(direct),
                 tolerance = 1e-10,
                 info = paste(row, collapse = " "))
  }
})

test_that("weighted non-categorical CV uses weighted row leverages", {
  data <- make_noncat_cv_gram_data(n = 180L, seed = 93001L)
  K <- matrix(c(3, 2, 2, 2), ncol = 2L, byrow = TRUE)

  for(basis in c("additive", "tensor", "glp")) {
    P <- crs:::prod.spline(
      x = data$x,
      K = K,
      knots = "quantiles",
      basis = basis,
      display.warnings = FALSE
    )
    X <- if(basis %in% c("additive", "glp")) cbind(1, P) else P

    fit <- crs:::.crs_weighted_ls_cv_rows(
      X = X,
      y = data$y,
      weights = data$weights,
      rows = seq_len(NROW(X)),
      ridge.lambda = NULL,
      rcond.min = 1e-10,
      allow.fallback = TRUE,
      use.svd.fallback = TRUE
    )

    correct_cv <- mean((sqrt(data$weights) * fit$residuals.rows)^2 /
                         (1 - fit$hat.rows)^2)
    legacy_cv <- mean((sqrt(data$weights) * fit$residuals.rows)^2 /
                        (1 - pmin(hat(P), 1 - .Machine$double.eps))^2)
    actual_cv <- noncat_cv_gram_value(data, basis, "cv.ls",
                                      weighted = TRUE, use_gram = TRUE)

    expect_equal(as.numeric(actual_cv), correct_cv, tolerance = 1e-12,
                 info = basis)
    expect_true(abs(legacy_cv - correct_cv) > 1e-8, info = basis)
  }
})
