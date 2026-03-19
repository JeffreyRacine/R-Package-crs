capture_messages_only <- function(expr) {
  messages <- character()
  withCallingHandlers(
    expr,
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  messages
}

normalize_messages <- function(x) {
  sub("\n$", "", x)
}

progress_time_values <- function(values) {
  force(values)
  i <- 0L
  function() {
    i <<- min(i + 1L, length(values))
    values[[i]]
  }
}

test_that("progress begin returns disabled state when messages are off", {
  begin <- getFromNamespace(".crs_progress_begin", "crs")

  old_opts <- options(crs.messages = FALSE)
  on.exit(options(old_opts), add = TRUE)

  state <- with_crs_progress_bindings(
    list(.crs_progress_is_interactive = function() TRUE,
         .crs_progress_now = function() 10),
    begin("NOMAD search")
  )

  expect_false(state$enabled)
  expect_identical(state$label, "NOMAD search")
  expect_null(state$total)
  expect_false(state$known_total)
})

test_that("unknown-total progress delays start note until grace elapses", {
  begin <- getFromNamespace(".crs_progress_begin", "crs")
  step <- getFromNamespace(".crs_progress_step", "crs")

  old_opts <- options(crs.messages = TRUE, crs.progress.start.grace.unknown.sec = 1.0)
  on.exit(options(old_opts), add = TRUE)

  messages <- with_crs_progress_bindings(
    list(
      .crs_progress_is_interactive = function() TRUE,
      .crs_progress_now = progress_time_values(c(0, 0.4, 1.2))
    ),
    capture_messages_only({
      state <- begin("NOMAD search", surface = "general")
      state <- step(state, done = 3, detail = "fv=1.2")
      step(state, done = 4, detail = "fv=0.9")
    })
  )

  messages <- normalize_messages(messages)

  expect_identical(messages[[1L]], "[crs] NOMAD search...")
  expect_identical(
    messages[[2L]],
    "[crs] NOMAD search... iteration 4, elapsed 1.2s: fv=0.9"
  )
  expect_length(messages, 2L)
})

test_that("single-line progress trace records NOMAD updates and finish", {
  begin <- getFromNamespace(".crs_progress_begin", "crs")
  step <- getFromNamespace(".crs_progress_step_at", "crs")
  finish <- getFromNamespace(".crs_progress_end", "crs")

  old_opts <- options(
    crs.messages = TRUE,
    crs.progress.start.grace.unknown.sec = 0,
    crs.progress.interval.unknown.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_crs_progress_shadow_trace(
    {
      state <- begin("NOMAD search", surface = "nomad")
      state <- step(state, now = 0, done = 1, detail = "fv=1.2", force = TRUE)
      state <- step(state, now = 2, done = 2, detail = "fv=0.8", force = TRUE)
      finish(state, detail = "best fv=0.8")
    },
    force_renderer = "single_line",
    now = progress_time_values(c(0, 0, 2, 2)),
    interactive = TRUE
  )

  expect_true(any(vapply(actual$trace, function(x) identical(x$event, "finish"), logical(1L))))
  expect_true(any(grepl("NOMAD search... iteration 2, elapsed 2.0s: fv=0.8", vapply(actual$trace, `[[`, character(1L), "line"), fixed = TRUE)))
})

test_that("plot progress emits checkpoint updates even when throttled", {
  begin <- getFromNamespace(".crs_plot_stage_progress_begin", "crs")
  tick <- getFromNamespace(".crs_plot_progress_tick", "crs")
  finish <- getFromNamespace(".crs_plot_progress_end", "crs")

  old_opts <- options(
    crs.messages = TRUE,
    crs.plot.progress.start.grace.sec = 0,
    crs.plot.progress.interval.sec = 999,
    crs.plot.progress.max.intermediate = 2
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_crs_progress_shadow_trace(
    {
      state <- begin(total = 8, label = "Plot bootstrap")
      state <- tick(state, done = 1)
      state <- tick(state, done = 3)
      state <- tick(state, done = 6)
      finish(state)
    },
    force_renderer = "single_line",
    now = progress_time_values(c(0, 0, 0.1, 0.2, 0.3, 0.4, 0.5)),
    interactive = TRUE
  )

  lines <- vapply(actual$trace, `[[`, character(1L), "line")
  expect_true(any(grepl("Plot bootstrap 0/8", lines, fixed = TRUE)))
  expect_true(any(grepl("Plot bootstrap 3/8", lines, fixed = TRUE)))
  expect_true(any(grepl("Plot bootstrap 6/8", lines, fixed = TRUE)))
})

test_that("console helpers route through the shared progress runtime", {
  reset <- getFromNamespace(".crs_progress_reset_registry", "crs")

  old_opts <- options(crs.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  with_crs_progress_bindings(
    list(
      .crs_progress_is_interactive = function() TRUE,
      .crs_progress_now = progress_time_values(c(0, 0, 0))
    ),
    {
      reset()
      console <- newLineConsole()
      console <- printPush("Working...", console)
      expect_false(is.null(console$progress$state))
      console <- printClear(console)
      expect_null(console$progress$state)
    }
  )
})
