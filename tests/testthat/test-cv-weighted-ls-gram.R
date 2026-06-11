make_cv_efficiency_data <- function(n = 90L, seed = 20260610L) {
  set.seed(seed)
  x1 <- runif(n)
  x2 <- runif(n)
  z_unordered <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  z_ordered <- ordered(sample(1:4, n, replace = TRUE))
  z_binary <- factor(sample(c("no", "yes"), n, replace = TRUE))
  y <- sin(2 * pi * x1) + 0.5 * cos(2 * pi * x2) +
    as.numeric(z_unordered) / 5 + as.numeric(z_ordered) / 7 +
    rnorm(n, sd = 0.15)
  list(
    x = cbind(x1 = x1, x2 = x2),
    z_factor = data.frame(
      z_unordered = z_unordered,
      z_ordered = z_ordered,
      z_binary = z_binary
    ),
    z_kernel = cbind(
      z_unordered = as.numeric(z_unordered),
      z_ordered = as.numeric(z_ordered),
      z_binary = as.numeric(z_binary)
    ),
    is_ordered = c(FALSE, TRUE, FALSE),
    y = y,
    weights = seq(0.6, 1.4, length.out = n)
  )
}

cv_efficiency_kernel_index <- function(z_kernel) {
  z_unique <- crs:::uniquecombs(z_kernel)
  ind <- attr(z_unique, "index")
  ind_vals <- unique(ind)
  ind_list <- lapply(ind_vals, function(i) ind == i)
  list(
    z_unique = z_unique,
    ind = ind,
    ind_vals = ind_vals,
    ind_list = ind_list,
    nrow_z_unique = NROW(z_unique)
  )
}

glp_reference_matrix <- function(X) {
  k <- length(X)
  dimen <- vapply(X, ncol, numeric(1L))
  dimen_list <- lapply(dimen, function(d) 0:d)
  z <- do.call("expand.grid", dimen_list)
  s <- rowSums(z)
  z <- z[(s > 0) & (s <= max(dimen)), , drop = FALSE]

  if(!all(dimen == max(dimen))) {
    for(j in seq_along(dimen)) {
      d <- dimen[[j]]
      if((d < max(dimen)) && (d > 0)) {
        s <- rowSums(z)
        drop <- (s > d) & (z[, j, drop = FALSE] ==
                             matrix(d, nrow(z), 1, byrow = TRUE))
        z <- z[!drop, , drop = FALSE]
      }
    }
  }

  out <- cbind(1, X[[1]])[, 1 + z[, 1]]
  if(k > 1L) {
    for(i in 2:k) {
      out <- out * cbind(1, X[[i]])[, 1 + z[, i]]
    }
  }
  matrix(out, nrow = NROW(X[[1]]))
}

qr_weighted_reference <- function(X, y, weights = NULL, ridge.lambda = NULL) {
  if(is.null(weights)) weights <- rep(1, NROW(X))
  sw <- sqrt(weights)
  if(!is.null(ridge.lambda)) {
    k <- NCOL(X)
    X_fit <- rbind(X, sqrt(ridge.lambda) * diag(k))
    y_fit <- c(y, rep(0, k))
    sw <- c(sw, rep(1, k))
  } else {
    X_fit <- X
    y_fit <- y
  }
  model <- .lm.fit(X_fit * sw, y_fit * sw, tol = 1e-7)
  fitted <- drop(X %*% model$coefficients)
  h <- crs:::hat.from.lm.fit(model)
  if(!is.null(ridge.lambda)) h <- h[seq_len(NROW(X))]
  list(coefficients = as.numeric(model$coefficients),
       residuals = y - fitted,
       hat = h)
}

run_kernel_cv_case <- function(data,
                               kidx,
                               basis,
                               cv_func = "cv.ls",
                               weights = NULL,
                               use_gram = TRUE,
                               gram_rcond_min = 1e-10) {
  K <- matrix(c(3, 2, 2, 2), ncol = 2L, byrow = TRUE)
  crs:::cv.kernel.spline(
    x = data$x,
    y = data$y,
    z = data$z_kernel,
    K = K,
    lambda = c(0.35, 0.6, 0.8),
    z.unique = kidx$z_unique,
    ind = kidx$ind,
    ind.vals = kidx$ind_vals,
    ind.list = kidx$ind_list,
    nrow.z.unique = kidx$nrow_z_unique,
    is.ordered.z = data$is_ordered,
    knots = "quantiles",
    basis = basis,
    cv.func = cv_func,
    weights = weights,
    tau = NULL,
    singular.ok = FALSE,
    display.warnings = FALSE,
    use.ridge = FALSE,
    smooth.penalty = TRUE,
    penalty.scale = 1000,
    use.gram.cv = use_gram,
    gram.rcond.min = gram_rcond_min,
    record.gram.stats = TRUE
  )
}

run_noncat_cv_case <- function(data,
                               basis,
                               cv_func = "cv.ls",
                               weights = NULL,
                               use_gram = TRUE,
                               gram_rcond_min = 1e-10) {
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
    weights = weights,
    tau = NULL,
    singular.ok = FALSE,
    display.warnings = FALSE,
    use.ridge = FALSE,
    smooth.penalty = TRUE,
    penalty.scale = 1000,
    use.gram.cv = use_gram,
    gram.rcond.min = gram_rcond_min,
    record.gram.stats = TRUE
  )
}

run_factor_cv_case <- function(data,
                               basis,
                               cv_func = "cv.ls",
                               weights = NULL,
                               use_gram = TRUE,
                               gram_rcond_min = 1e-10) {
  K <- matrix(c(2, 1, 2, 1), ncol = 2L, byrow = TRUE)
  crs:::cv.factor.spline(
    x = data$x,
    y = data$y,
    z = data$z_factor,
    K = K,
    I = c(1, 1, 1),
    knots = "quantiles",
    basis = basis,
    cv.func = cv_func,
    weights = weights,
    tau = NULL,
    singular.ok = FALSE,
    display.warnings = FALSE,
    use.ridge = FALSE,
    smooth.penalty = TRUE,
    penalty.scale = 1000,
    use.gram.cv = use_gram,
    gram.rcond.min = gram_rcond_min,
    record.gram.stats = TRUE
  )
}

test_that("glp.model.matrix preserves the filtered GLP column oracle", {
  set.seed(20260610)
  dims_cases <- list(c(1, 1), c(2, 1), c(2, 2), c(3, 2, 1), c(4, 2, 3))

  for(dims in dims_cases) {
    X <- lapply(dims, function(d) matrix(rnorm(45L * d), ncol = d))
    actual <- crs:::glp.model.matrix(X)
    expected <- glp_reference_matrix(X)
    expect_equal(dim(actual), dim(expected))
    expect_equal(actual, expected, tolerance = 0)
    expect_equal(qr(actual)$rank, qr(expected)$rank)
  }
})

test_that("weighted LS primitive matches QR residual and leverage references", {
  data <- make_cv_efficiency_data(n = 110L, seed = 93001L)
  K <- matrix(c(3, 2, 2, 2), ncol = 2L, byrow = TRUE)
  rows <- sort(unique(c(1L, seq(3L, NROW(data$x), by = 9L), NROW(data$x))))

  cases <- expand.grid(
    basis = c("additive", "tensor", "glp"),
    weighted = c(FALSE, TRUE),
    ridge = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )

  for(i in seq_len(nrow(cases))) {
    row <- cases[i, ]
    P <- crs:::prod.spline(data$x, K = K, knots = "quantiles",
                           basis = row$basis, display.warnings = FALSE)
    X <- if(row$basis %in% c("additive", "glp")) cbind(1, P) else P
    weights <- if(row$weighted) data$weights else NULL
    ridge.lambda <- if(row$ridge) 1e-3 else NULL
    ref <- qr_weighted_reference(X, data$y, weights, ridge.lambda)

    core <- crs:::.crs_weighted_ls_core(
      X, data$y, weights = weights, ridge.lambda = ridge.lambda,
      rcond.min = 0, allow.fallback = FALSE
    )
    expect_identical(core$status, "gram")
    expect_equal(as.numeric(core$coefficients), ref$coefficients,
                 tolerance = 1e-8)

    cv_rows <- crs:::.crs_weighted_ls_cv_rows(
      X, data$y, weights = weights, rows = rows,
      ridge.lambda = ridge.lambda, rcond.min = 0,
      allow.fallback = FALSE
    )
    expect_equal(cv_rows$residuals.rows, ref$residuals[rows],
                 tolerance = 1e-8)
    expect_equal(cv_rows$hat.rows, ref$hat[rows], tolerance = 1e-8)

    fallback_rows <- crs:::.crs_weighted_ls_cv_rows(
      X, data$y, weights = weights, rows = rows,
      ridge.lambda = ridge.lambda, rcond.min = 1,
      allow.fallback = TRUE
    )
    expect_identical(fallback_rows$status, "fallback_rcond")
    expect_true(fallback_rows$method %in% c("qr", "svd"))
    expect_equal(fallback_rows$residuals.rows, ref$residuals[rows],
                 tolerance = 1e-8)
    expect_equal(fallback_rows$hat.rows, ref$hat[rows], tolerance = 1e-8)
  }
})

test_that("QR fallback residuals are pivot-safe for rank-deficient designs", {
  set.seed(20260611)
  n <- 18L
  x <- seq_len(n) / n
  X <- cbind(
    intercept = 1,
    x = x,
    x_duplicate = x,
    z = rep(c(0, 1), length.out = n)
  )
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.05)
  rows <- c(1L, 4L, 9L, 13L, 18L)

  cases <- expand.grid(
    weighted = c(FALSE, TRUE),
    ridge = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )

  for(i in seq_len(nrow(cases))) {
    row <- cases[i, ]
    weights <- if(row$weighted) seq(0.6, 1.4, length.out = n) else NULL
    ridge.lambda <- if(row$ridge) 1e-3 else NULL
    weights.fit <- if(is.null(weights)) rep(1, n) else weights
    sw <- sqrt(weights.fit)

    if(is.null(ridge.lambda)) {
      X.fit <- X
      y.fit <- y
      sw.fit <- sw
    } else {
      p <- ncol(X)
      X.fit <- rbind(X, sqrt(ridge.lambda) * diag(p))
      y.fit <- c(y, rep(0, p))
      sw.fit <- c(sw, rep(1, p))
    }

    ref <- .lm.fit(X.fit * sw.fit, y.fit * sw.fit, tol = 1e-7)
    ref.residuals <- ref$residuals[seq_len(n)] / sw
    ref.hat <- crs:::hat.from.lm.fit(ref)
    ref.hat <- ref.hat[seq_len(n)]

    fit <- crs:::.crs_weighted_ls_cv_rows(
      X = X,
      y = y,
      weights = weights,
      rows = rows,
      ridge.lambda = ridge.lambda,
      rcond.min = 1,
      allow.fallback = TRUE
    )

    expect_true(fit$status %in% c("fallback_rcond", "fallback_chol"))
    expect_true(fit$method %in% c("qr", "svd"))
    expect_equal(fit$residuals.rows, ref.residuals[rows], tolerance = 1e-12)
    expect_equal(fit$hat.rows, ref.hat[rows], tolerance = 1e-12)
    expect_false(anyNA(fit$residuals.rows))
  }
})

test_that("categorical kernel CV Gram path preserves objectives and fallback", {
  data <- make_cv_efficiency_data(n = 120L, seed = 94001L)
  kidx <- cv_efficiency_kernel_index(data$z_kernel)
  cases <- expand.grid(
    basis = c("additive", "tensor", "glp"),
    cv_func = c("cv.ls", "cv.gcv", "cv.aic"),
    weighted = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )

  for(i in seq_len(nrow(cases))) {
    row <- cases[i, ]
    weights <- if(row$weighted) data$weights else NULL
    current <- run_kernel_cv_case(data, kidx, row$basis, row$cv_func,
                                  weights, use_gram = FALSE)
    gram <- run_kernel_cv_case(data, kidx, row$basis, row$cv_func,
                               weights, use_gram = TRUE)
    fallback <- run_kernel_cv_case(data, kidx, row$basis, row$cv_func,
                                   weights, use_gram = TRUE,
                                   gram_rcond_min = 1)
    gram_stats <- attr(gram, "gram.stats", exact = TRUE)
    fallback_stats <- attr(fallback, "gram.stats", exact = TRUE)
    expect_equal(as.numeric(gram), as.numeric(current), tolerance = 1e-8)
    expect_equal(as.numeric(fallback), as.numeric(current), tolerance = 1e-8)
    expect_gt(gram_stats$used + gram_stats$fallback_rcond +
                gram_stats$fallback_chol, 0)
    expect_equal(fallback_stats$used, 0)
    expect_gt(fallback_stats$fallback_rcond, 0)
  }
})

test_that("non-categorical kernel CV Gram path preserves objectives and fallback", {
  data <- make_cv_efficiency_data(n = 120L, seed = 95001L)
  cases <- expand.grid(
    basis = c("additive", "tensor", "glp"),
    cv_func = c("cv.ls", "cv.gcv", "cv.aic"),
    weighted = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )

  for(i in seq_len(nrow(cases))) {
    row <- cases[i, ]
    weights <- if(row$weighted) data$weights else NULL
    current <- run_noncat_cv_case(data, row$basis, row$cv_func, weights,
                                  use_gram = FALSE)
    gram <- run_noncat_cv_case(data, row$basis, row$cv_func, weights,
                               use_gram = TRUE)
    fallback <- run_noncat_cv_case(data, row$basis, row$cv_func, weights,
                                   use_gram = TRUE, gram_rcond_min = 1)
    gram_stats <- attr(gram, "gram.stats", exact = TRUE)
    fallback_stats <- attr(fallback, "gram.stats", exact = TRUE)
    expect_equal(as.numeric(gram), as.numeric(current), tolerance = 1e-8)
    expect_equal(as.numeric(fallback), as.numeric(current), tolerance = 1e-8)
    expect_gt(gram_stats$used + gram_stats$fallback_rcond +
                gram_stats$fallback_chol, 0)
    expect_equal(fallback_stats$used, 0)
    expect_gt(fallback_stats$fallback_rcond, 0)
  }
})

test_that("factor CV Gram path preserves objectives and fallback", {
  data <- make_cv_efficiency_data(n = 180L, seed = 96001L)
  cases <- expand.grid(
    basis = c("additive", "tensor", "glp"),
    cv_func = c("cv.ls", "cv.gcv", "cv.aic"),
    weighted = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )

  for(i in seq_len(nrow(cases))) {
    row <- cases[i, ]
    weights <- if(row$weighted) data$weights else NULL
    current <- run_factor_cv_case(data, row$basis, row$cv_func, weights,
                                  use_gram = FALSE)
    gram <- run_factor_cv_case(data, row$basis, row$cv_func, weights,
                               use_gram = TRUE)
    fallback <- run_factor_cv_case(data, row$basis, row$cv_func, weights,
                                   use_gram = TRUE, gram_rcond_min = 1)
    gram_stats <- attr(gram, "gram.stats", exact = TRUE)
    fallback_stats <- attr(fallback, "gram.stats", exact = TRUE)
    expect_equal(as.numeric(gram), as.numeric(current), tolerance = 1e-8)
    expect_equal(as.numeric(fallback), as.numeric(current), tolerance = 1e-8)
    expect_gt(gram_stats$used + gram_stats$fallback_rcond +
                gram_stats$fallback_chol, 0)
    expect_equal(fallback_stats$used, 0)
    expect_gt(fallback_stats$fallback_rcond, 0)
  }
})
