test_that("snomadr rejects unsupported NOMAD OpenMP-only options", {
  eval_f <- function(x) sum((x - 0.25)^2)

  coop_messages <- capture.output(
    coop_result <- crs::snomadr(
      n = 2,
      eval.f = eval_f,
      x0 = c(0, 0),
      lb = c(-1, -1),
      ub = c(1, 1),
      opts = list("COOP_MADS_OPTIMIZATION" = "yes"),
      display.nomad.progress = FALSE
    ),
    type = "message"
  )
  expect_equal(coop_result$status, 1)
  expect_match(coop_result$message, "invalid parameters")
  expect_true(any(grepl(
    "COOP_MADS_OPTIMIZATION is not supported in the shipped crs build",
    coop_messages,
    fixed = TRUE
  )))

  psd_messages <- capture.output(
    psd_result <- crs::snomadr(
      n = 2,
      eval.f = eval_f,
      x0 = c(0, 0),
      lb = c(-1, -1),
      ub = c(1, 1),
      opts = list("PSD_MADS_OPTIMIZATION" = "yes"),
      display.nomad.progress = FALSE
    ),
    type = "message"
  )
  expect_equal(psd_result$status, 1)
  expect_match(psd_result$message, "invalid parameters")
  expect_true(any(grepl(
    "PSD_MADS_OPTIMIZATION is not supported in the shipped crs build",
    psd_messages,
    fixed = TRUE
  )))

  thread_messages <- capture.output(
    thread_result <- crs::snomadr(
      n = 2,
      eval.f = eval_f,
      x0 = c(0, 0),
      lb = c(-1, -1),
      ub = c(1, 1),
      opts = list("NB_THREADS_PARALLEL_EVAL" = "2"),
      display.nomad.progress = FALSE
    ),
    type = "message"
  )
  expect_equal(thread_result$status, 1)
  expect_match(thread_result$message, "invalid parameters")
  expect_true(any(grepl(
    "NB_THREADS_PARALLEL_EVAL>1 is not supported in the shipped crs build",
    thread_messages,
    fixed = TRUE
  )))
})

test_that("snomadr still succeeds on the shipped serial NOMAD path", {
  eval_f <- function(x) sum((x - 0.25)^2)

  result <- crs::snomadr(
    n = 2,
    eval.f = eval_f,
    x0 = c(0, 0),
    lb = c(-1, -1),
    ub = c(1, 1),
    opts = list(
      "MAX_BB_EVAL" = "15",
      "NB_THREADS_PARALLEL_EVAL" = "1"
    ),
    display.nomad.progress = FALSE
  )

  expect_type(result, "list")
  expect_true(is.numeric(result$objective))
  expect_length(result$solution, 2)
})
