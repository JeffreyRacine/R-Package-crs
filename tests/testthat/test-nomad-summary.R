test_that("NOMAD summary prints compact native cache diagnostics", {
  print_summary <- getFromNamespace(".crs_nomad_summary_print", "crs")
  object <- list(nomad.summary = list(
    num.feval = 252,
    cache.hits = 19762,
    total.evaluations = 20014,
    callback.evaluations = 252,
    cache.size = 252,
    callback.cache.hits = NA_real_,
    callback.cache.misses = NA_real_,
    callback.cache.size = NA_real_
  ))

  out <- capture.output(print_summary(object))

  expect_true(any(grepl("Number of Function Evaluations: 252",
                        out, fixed = TRUE)))
  expect_true(any(grepl("NOMAD cache: 98.7% hits (19,762/20,014 point requests)",
                        out, fixed = TRUE)))
  expect_false(any(grepl("total = ", out, fixed = TRUE)))
  expect_false(any(grepl("callbacks = ", out, fixed = TRUE)))
  expect_false(any(grepl("NOMAD Cache Summary", out, fixed = TRUE)))
})

test_that("NOMAD summary omits native cache line when there are no hits", {
  print_summary <- getFromNamespace(".crs_nomad_summary_print", "crs")
  object <- list(nomad.summary = list(
    num.feval = 316,
    cache.hits = 0,
    total.evaluations = 316,
    callback.evaluations = 316,
    callback.cache.hits = NA_real_,
    callback.cache.misses = NA_real_,
    callback.cache.size = NA_real_
  ))

  out <- capture.output(print_summary(object))

  expect_true(any(grepl("Number of Function Evaluations: 316",
                        out, fixed = TRUE)))
  expect_false(any(grepl("NOMAD cache:", out, fixed = TRUE)))
})

test_that("NOMAD summary handles partial native cache metadata conservatively", {
  print_summary <- getFromNamespace(".crs_nomad_summary_print", "crs")
  object <- list(nomad.summary = list(
    num.feval = 10,
    cache.hits = 1200,
    total.evaluations = NA_real_,
    callback.evaluations = 10,
    callback.cache.hits = NA_real_,
    callback.cache.misses = NA_real_,
    callback.cache.size = NA_real_
  ))

  out <- capture.output(print_summary(object))

  expect_true(any(grepl("Number of Function Evaluations: 10",
                        out, fixed = TRUE)))
  expect_true(any(grepl("NOMAD cache: 1,200 hits",
                        out, fixed = TRUE)))
  expect_false(any(grepl("point requests", out, fixed = TRUE)))
})

test_that("NOMAD summary prints CRS restart cache separately", {
  print_summary <- getFromNamespace(".crs_nomad_summary_print", "crs")
  object <- list(nomad.summary = list(
    num.feval = 252,
    cache.hits = 19762,
    total.evaluations = 20014,
    callback.evaluations = 252,
    callback.cache.hits = 14,
    callback.cache.misses = 238,
    callback.cache.size = 238
  ))

  out <- capture.output(print_summary(object))

  expect_true(any(grepl("NOMAD cache: 98.7% hits (19,762/20,014 point requests)",
                        out, fixed = TRUE)))
  expect_true(any(grepl("CRS restart cache: 14 hits, 238 misses, 238 stored points",
                        out, fixed = TRUE)))
  expect_false(any(grepl("NOMAD Cache Summary", out, fixed = TRUE)))
})

test_that("NOMAD summary printer is used by live CRS summaries", {
  set.seed(6401)
  n <- 36L
  x <- runif(n)
  z <- factor(sample(c("a", "b"), n, replace = TRUE))
  y <- sin(2 * pi * x) + 0.2 * (z == "b") + rnorm(n, sd = 0.04)

  fit <- suppressWarnings(crs(
    y ~ x + z,
    kernel = TRUE,
    basis = "additive",
    cv = "nomad",
    degree.max = 2,
    segments.max = 2,
    nmulti = 1,
    max.bb.eval = 12,
    random.seed = 6401,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  ))

  out <- capture.output(summary(fit))

  expect_true(any(grepl("Number of Function Evaluations:", out, fixed = TRUE)))
  expect_false(any(grepl("total = ", out, fixed = TRUE)))
  expect_false(any(grepl("callbacks = ", out, fixed = TRUE)))
})
