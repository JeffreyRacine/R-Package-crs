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

capture_single_line_output <- function(bindings, code) {
  path <- tempfile(fileext = ".log")
  con <- file(path, open = "wb")
  closed <- FALSE
  reset <- getFromNamespace(".crs_progress_reset_registry", "crs")

  on.exit({
    if (!closed) {
      close(con)
    }
    unlink(path)
  }, add = TRUE)

  bindings$.crs_progress_single_line_connection <- function() con
  reset()
  on.exit(reset(), add = TRUE)
  with_crs_progress_bindings(bindings, code)

  close(con)
  closed <- TRUE
  readChar(path, nchars = file.info(path)$size, useBytes = TRUE)
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

test_that("status helpers route through the shared progress runtime", {
  reset <- getFromNamespace(".crs_progress_reset_registry", "crs")
  begin_status <- getFromNamespace(".crs_progress_status_begin", "crs")
  update_status <- getFromNamespace(".crs_progress_status_update", "crs")
  clear_status <- getFromNamespace(".crs_progress_status_clear", "crs")

  old_opts <- options(crs.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  with_crs_progress_bindings(
    list(
      .crs_progress_is_interactive = function() TRUE,
      .crs_progress_now = progress_time_values(c(0, 0, 0))
    ),
    {
      reset()
      status <- begin_status(surface = "console")
      status <- update_status(status, "Working...")
      expect_false(is.null(status$state))
      status <- clear_status(status)
      expect_null(status$state)
    }
  )
})

test_that("manual status updates honor the shared throttle and retain the latest line", {
  begin_status <- getFromNamespace(".crs_progress_status_begin", "crs")
  update_status <- getFromNamespace(".crs_progress_status_update", "crs")
  clear_status <- getFromNamespace(".crs_progress_status_clear", "crs")

  old_opts <- options(
    crs.messages = TRUE,
    crs.progress.interval.unknown.sec = 2
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_crs_progress_shadow_trace(
    {
      status <- begin_status(surface = "nomad")
      update_status(status, "Selecting spline model... iteration 1, elapsed 0.0s: fv=1.0")
      update_status(status, "Selecting spline model... iteration 2, elapsed 0.2s: fv=0.8")
      update_status(status, "Selecting spline model... iteration 3, elapsed 0.5s: fv=0.6")
      clear_status(status)
    },
    force_renderer = "single_line",
    now = progress_time_values(c(0, 0, 0.2, 0.5, 0.5)),
    interactive = TRUE
  )

  lines <- vapply(actual$trace, `[[`, character(1L), "line")
  events <- vapply(actual$trace, `[[`, character(1L), "event")

  expect_identical(sum(events == "render"), 1L)
  expect_true(any(grepl("iteration 1", lines, fixed = TRUE)))
  expect_true(any(grepl("iteration 3", lines, fixed = TRUE)))
  expect_true(identical(events[[length(events)]], "finish"))
})

test_that("manual status intervals can differ by surface", {
  begin_status <- getFromNamespace(".crs_progress_status_begin", "crs")
  update_status <- getFromNamespace(".crs_progress_status_update", "crs")

  old_opts <- options(
    crs.messages = TRUE,
    crs.progress.interval.unknown.sec = 5,
    crs.progress.manual.interval.nomad.sec = 5,
    crs.progress.manual.interval.cv.sec = 1
  )
  on.exit(options(old_opts), add = TRUE)

  nomad <- begin_status(surface = "nomad")
  cv <- begin_status(surface = "cv")

  update_status(nomad, "nomad iteration 1")
  update_status(cv, "cv iteration 1")

  expect_identical(nomad$state$throttle_sec, 5)
  expect_identical(cv$state$throttle_sec, 1)
})

test_that("IV wrappers can suppress pending NOMAD start frames without losing first real eval", {
  begin_status <- getFromNamespace(".crs_progress_status_begin", "crs")
  clear_status <- getFromNamespace(".crs_progress_status_clear", "crs")
  begin_nomad <- getFromNamespace(".crs_nomad_progress_begin", "crs")
  maybe_restart <- getFromNamespace(".crs_nomad_progress_maybe_restart", "crs")
  eval_nomad <- getFromNamespace(".crs_nomad_progress_eval", "crs")

  old_opts <- options(
    crs.messages = TRUE,
    crs.progress.interval.unknown.sec = 5,
    crs.nomad.progress.show.start = FALSE
  )
  on.exit(options(old_opts), add = TRUE)

  traced <- capture_crs_progress_shadow_trace(
    {
      status <- begin_status(surface = "nomad")
      progress <- begin_nomad(
        progress.status = status,
        label = "Selecting spline model",
        starts = matrix(c(1, 1), nrow = 1L)
      )
      maybe_restart(
        progress = progress,
        input = c(1, 1),
        degree = 1L,
        segments = 1L
      )
      eval_nomad(
        progress = progress,
        degree = 1L,
        segments = 1L,
        value = 0.25
      )
      clear_status(status)
    },
    force_renderer = "legacy",
    now = progress_time_values(c(0, 0, 0.1, 0.2)),
    interactive = TRUE
  )

  lines <- unique(vapply(traced$trace, `[[`, character(1L), "line"))
  expect_false(any(grepl("fv=pending", lines, fixed = TRUE)))
  expect_false(any(grepl("best deg pending", lines, fixed = TRUE)))
  expect_true(any(grepl("iteration 1", lines, fixed = TRUE)))
  expect_true(any(grepl("eval 1", lines, fixed = TRUE)))
  expect_true(any(grepl("fv=0.25", lines, fixed = TRUE)))
})

test_that("single-line NOMAD rendering keeps detail when the line is wide", {
  begin_status <- getFromNamespace(".crs_progress_status_begin", "crs")
  update_status <- getFromNamespace(".crs_progress_status_update", "crs")
  clear_status <- getFromNamespace(".crs_progress_status_clear", "crs")
  rendered <- character()

  old_opts <- options(
    crs.messages = TRUE,
    crs.progress.interval.unknown.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  with_crs_progress_bindings(
    list(
      .crs_progress_output_width = function() 70L,
      .crs_progress_is_interactive = function() TRUE,
      .crs_progress_now = progress_time_values(c(0, 0, 0)),
      .crs_progress_render_single_line = function(snapshot, event = c("render", "finish", "abort")) {
        rendered <<- c(rendered, snapshot$render_line)
        invisible(snapshot)
      }
    ),
    {
      status <- begin_status(surface = "nomad")
      update_status(
        status,
        "Selecting spline model... iteration 10, elapsed 5.0s: multistart 2/2, eval 16, deg (5,6), seg (9,2), best deg (4,7), best seg (10,1), include (1,0), fv=0.5008949"
      )
      clear_status(status)
    }
  )

  expect_true(any(grepl("include", rendered, fixed = TRUE)))
  expect_true(any(grepl("fv=", rendered, fixed = TRUE)))
  expect_false(any(grepl("^\\[crs\\] Selecting spline model\\.\\.\\. iteration 10, elapsed 5\\.0s$", rendered)))
})

test_that("single-line finish clear does not emit a trailing newline", {
  render <- getFromNamespace(".crs_progress_render_single_line", "crs")

  output <- capture_single_line_output(
    list(.crs_progress_output_width = function() 24L),
    {
      render(
        list(
          render_line = "[crs] Working...",
          last_width = 18L
        ),
        event = "finish"
      )
    }
  )

  expect_false(grepl("\n", output, fixed = TRUE))
  expect_true(startsWith(output, "\r"))
  expect_true(endsWith(output, "\r"))
})

test_that("single-line finish clear uses ANSI erase on capable ttys", {
  render <- getFromNamespace(".crs_progress_render_single_line", "crs")

  output <- capture_single_line_output(
    list(.crs_progress_single_line_supports_ansi = function(con) TRUE),
    {
      render(
        list(
          render_line = "[crs] Working...",
          last_width = 18L
        ),
        event = "finish"
      )
    }
  )

  expect_false(grepl("\n", output, fixed = TRUE))
  expect_identical(output, "\r\033[2K\r")
})

test_that("single-line render clears stale suffix across state handoff in non-ANSI mode", {
  render <- getFromNamespace(".crs_progress_render_single_line", "crs")
  output_width <- 120L
  long <- "[crs] Calling NOMAD (Nonsmooth Optimization by Mesh Adaptive Direct Search)... elapsed 0.0s"
  short <- "[crs] NOMAD search..."
  medium <- "[crs] NOMAD search... iteration 20, elapsed 2.5s: fv=0.2585065"

  output <- capture_single_line_output(
    list(
      .crs_progress_output_width = function() output_width,
      .crs_progress_single_line_supports_ansi = function(con) FALSE
    ),
    {
      render(list(render_line = long, last_width = 0L), event = "render")
      render(list(render_line = short, last_width = 0L), event = "render")
      render(list(render_line = medium, last_width = nchar(short, type = "width")), event = "render")
    }
  )

  expect_identical(
    output,
    paste0(
      "\r", long,
      "\r", short, strrep(" ", nchar(long, type = "width") - nchar(short, type = "width")),
      "\r", medium, strrep(" ", nchar(long, type = "width") - nchar(medium, type = "width"))
    )
  )
})

test_that("activity helpers emit and clear through the shared renderer", {
  begin_activity <- getFromNamespace(".crs_progress_activity_begin", "crs")
  end_activity <- getFromNamespace(".crs_progress_activity_end", "crs")

  old_opts <- options(crs.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_crs_progress_shadow_trace(
    {
      state <- begin_activity("Calling NOMAD", surface = "nomad")
      end_activity(state)
    },
    force_renderer = "single_line",
    now = progress_time_values(c(0, 0, 0)),
    interactive = TRUE
  )

  expect_true(any(vapply(actual$trace, function(x) identical(x$event, "finish"), logical(1L))))
  expect_true(any(grepl("Calling NOMAD", vapply(actual$trace, `[[`, character(1L), "line"), fixed = TRUE)))
})
