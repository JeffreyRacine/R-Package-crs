make_cell_cache_data <- function(kind = c("balanced", "unbalanced", "sparse"),
                                 n = 180L,
                                 seed = 20260611L) {
  kind <- match.arg(kind)
  set.seed(seed)
  x1 <- runif(n)
  x2 <- runif(n)

  if(kind == "balanced") {
    z1 <- factor(sample(letters[1:3], n, replace = TRUE))
    z2 <- ordered(sample(1:3, n, replace = TRUE))
  } else if(kind == "unbalanced") {
    z1 <- factor(sample(letters[1:4], n, replace = TRUE,
                        prob = c(0.55, 0.25, 0.15, 0.05)))
    z2 <- ordered(sample(1:4, n, replace = TRUE,
                         prob = c(0.55, 0.25, 0.15, 0.05)))
  } else {
    z1 <- factor(sample(sprintf("u%02d", 1:10), n, replace = TRUE))
    z2 <- ordered(sample(1:10, n, replace = TRUE))
  }

  z <- cbind(z1 = as.numeric(z1), z2 = as.numeric(z2))
  signal <- sin(2 * pi * x1) + 0.6 * cos(2 * pi * x2) +
    as.numeric(z1) / max(as.numeric(z1)) +
    0.4 * as.numeric(z2) / max(as.numeric(z2))
  list(
    x = cbind(x1 = x1, x2 = x2),
    y = signal + rnorm(n, sd = 0.2),
    z = z,
    is_ordered = c(FALSE, TRUE),
    weights = seq(0.6, 1.4, length.out = n)
  )
}

cell_cache_kernel_index <- function(z) {
  z_unique <- crs:::uniquecombs(z)
  ind <- attr(z_unique, "index")
  ind_vals <- unique(ind)
  list(
    z_unique = z_unique,
    ind = ind,
    ind_vals = ind_vals,
    ind_list = lapply(ind_vals, function(i) ind == i),
    nrow_z_unique = NROW(z_unique)
  )
}

cell_cache_cv <- function(data,
                          kidx,
                          basis,
                          cv_func,
                          weighted,
                          use_cell_cache,
                          use_gram = TRUE) {
  K <- matrix(c(3, 2, 2, 2), ncol = 2L, byrow = TRUE)
  crs:::cv.kernel.spline(
    x = data$x,
    y = data$y,
    z = data$z,
    K = K,
    lambda = c(0.4, 0.35),
    z.unique = kidx$z_unique,
    ind = kidx$ind,
    ind.vals = kidx$ind_vals,
    ind.list = kidx$ind_list,
    nrow.z.unique = kidx$nrow_z_unique,
    is.ordered.z = data$is_ordered,
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
    record.gram.stats = TRUE,
    use.cell.cache = use_cell_cache
  )
}

test_that("categorical kernel cell cache preserves CV objectives", {
  cases <- expand.grid(
    kind = c("balanced", "unbalanced", "sparse"),
    basis = c("additive", "tensor", "glp"),
    cv_func = c("cv.ls", "cv.gcv", "cv.aic"),
    weighted = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )

  for(i in seq_len(nrow(cases))) {
    row <- cases[i, ]
    data <- make_cell_cache_data(row$kind, seed = 71000L + i)
    kidx <- cell_cache_kernel_index(data$z)

    direct <- cell_cache_cv(data, kidx, row$basis, row$cv_func,
                            row$weighted, use_cell_cache = FALSE)
    cached <- cell_cache_cv(data, kidx, row$basis, row$cv_func,
                            row$weighted, use_cell_cache = TRUE)

    expect_equal(as.numeric(cached), as.numeric(direct),
                 tolerance = 1e-10,
                 info = paste(row, collapse = " "))

    cached_stats <- attr(cached, "gram.stats")
    expect_equal(cached_stats$fallback_chol, 0L)
    expect_equal(cached_stats$fallback_rcond, 0L)
    expect_equal(cached_stats$used, length(kidx$ind_vals))
  }
})

test_that("cell cache toggle does not override disabled Gram CV route", {
  data <- make_cell_cache_data("balanced", seed = 72001L)
  kidx <- cell_cache_kernel_index(data$z)

  direct <- cell_cache_cv(data, kidx, "additive", "cv.ls", FALSE,
                          use_cell_cache = FALSE, use_gram = FALSE)
  cached_arg <- cell_cache_cv(data, kidx, "additive", "cv.ls", FALSE,
                              use_cell_cache = TRUE, use_gram = FALSE)

  expect_equal(as.numeric(cached_arg), as.numeric(direct), tolerance = 0)
  expect_null(attr(cached_arg, "gram.stats"))
})

test_that("batched cell cache systems match per-cell recombination", {
  data <- make_cell_cache_data("sparse", seed = 73001L)
  kidx <- cell_cache_kernel_index(data$z)
  P <- crs:::prod.spline(
    x = data$x,
    K = matrix(c(3, 2, 2, 2), ncol = 2L, byrow = TRUE),
    knots = "quantiles",
    basis = "glp",
    display.warnings = FALSE
  )
  X <- cbind(1, P)
  cache <- crs:::.crs_cell_cache_build(
    X = X,
    y = data$y,
    ind = kidx$ind,
    ind.vals = kidx$ind_vals,
    weights = data$weights
  )
  cell_weights <- crs:::.crs_kernel_cell_weights(
    z.unique = kidx$z_unique,
    ind.vals = kidx$ind_vals,
    lambda = c(0.35, 0.65),
    is.ordered.z = data$is_ordered
  )
  systems <- crs:::.crs_cell_cache_systems(cache, cell_weights, max.bytes = Inf)

  expect_false(is.null(systems))
  for(i in seq_along(kidx$ind_vals)) {
    direct <- crs:::.crs_cell_cache_system(cache, cell_weights[i, ])
    batched <- crs:::.crs_cell_cache_system_at(systems, i)
    expect_equal(batched$G, direct$G, tolerance = 1e-12)
    expect_equal(batched$rhs, direct$rhs, tolerance = 1e-12)
  }

  expect_null(crs:::.crs_cell_cache_systems(cache, cell_weights, max.bytes = 1))
})
