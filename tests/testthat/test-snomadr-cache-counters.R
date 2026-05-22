test_that("snomadr exposes NOMAD native cache counters", {
  run_case <- function(eval_use_cache) {
    env <- new.env(parent = baseenv())
    env$count <- 0L
    env$target <- c(2, 3)

    eval.f <- function(x) {
      count <<- count + 1L
      sum((x - target)^2)
    }

    x0 <- matrix(
      c(1, 1,
        1, 1,
        1, 1,
        1, 1),
      nrow = 4,
      byrow = TRUE
    )

    out <- utils::capture.output({
      fit <- crs::snomadr(
        eval.f = eval.f,
        n = 2L,
        x0 = x0,
        bbin = c(1L, 1L),
        lb = c(1, 1),
        ub = c(5, 5),
        nmulti = 4L,
        random.seed = 42L,
        opts = list(
          DISPLAY_DEGREE = 0,
          MAX_BB_EVAL = 20,
          EVAL_USE_CACHE = eval_use_cache
        ),
        display.nomad.progress = FALSE,
        snomadr.environment = env
      )
    })
    invisible(out)

    list(fit = fit, callback.count = env$count)
  }

  cache_on <- run_case(TRUE)
  cache_off <- run_case(FALSE)

  expect_named(
    cache_on$fit,
    c("eval.f", "n", "bbin", "bbout", "x0", "lower.bounds",
      "upper.bounds", "nmulti", "random.seed", "options", "print.output",
      "snomadr.environment", "call", "status", "message", "bbe",
      "cache.hits", "cache.size", "callback.evaluations",
      "total.evaluations", "iterations", "objective", "solution")
  )
  expect_gt(cache_on$fit$cache.hits, 0L)
  expect_gte(cache_on$fit$cache.size, cache_on$fit$callback.evaluations)
  expect_equal(cache_on$fit$callback.evaluations, cache_on$callback.count)
  expect_equal(
    cache_on$fit$total.evaluations,
    cache_on$fit$callback.evaluations + cache_on$fit$cache.hits
  )

  expect_equal(cache_off$fit$cache.hits, 0L)
  expect_equal(cache_off$fit$callback.evaluations, cache_off$callback.count)
  expect_equal(
    cache_off$fit$total.evaluations,
    cache_off$fit$callback.evaluations
  )
})
