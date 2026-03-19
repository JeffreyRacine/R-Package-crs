.crs_progress_is_interactive <- function() {
  interactive()
}

.crs_progress_now <- function() {
  as.numeric(proc.time()[["elapsed"]])
}

.crs_progress_enabled <- function(domain = "general") {
  isTRUE(getOption("crs.messages", TRUE)) &&
    isTRUE(.crs_progress_is_interactive())
}

.crs_io_pkg_prefix <- function() {
  "[crs]"
}

.crs_io_prefix_text <- function(text) {
  prefix <- .crs_io_pkg_prefix()
  if (!is.character(text)) {
    return(text)
  }

  vapply(
    text,
    function(line) {
      if (is.na(line) || startsWith(line, prefix)) {
        return(line)
      }

      paste(prefix, line)
    },
    character(1L)
  )
}

.crs_message <- function(...) {
  message(.crs_io_prefix_text(paste(..., sep = "", collapse = "")))
  invisible(NULL)
}

.crs_warning <- function(..., call. = TRUE, immediate. = FALSE, noBreaks. = FALSE) {
  warning(
    .crs_io_prefix_text(paste(..., sep = "", collapse = "")),
    call. = call.,
    immediate. = immediate.,
    noBreaks. = noBreaks.
  )
}

.crs_progress_call_name <- function(call) {
  head <- call[[1L]]

  if (is.symbol(head)) {
    return(as.character(head))
  }

  if (is.call(head) &&
      identical(head[[1L]], as.name("::")) &&
      length(head) >= 3L &&
      is.symbol(head[[2L]]) &&
      is.symbol(head[[3L]])) {
    return(paste(as.character(head[[2L]]), as.character(head[[3L]]), sep = "::"))
  }

  ""
}

.crs_progress_is_message_muffled <- function() {
  calls <- sys.calls()
  if (!length(calls)) {
    return(FALSE)
  }

  any(vapply(
    calls,
    function(call) {
      .crs_progress_call_name(call) %in% c(
        "suppressMessages",
        "base::suppressMessages",
        "suppressPackageStartupMessages",
        "base::suppressPackageStartupMessages"
      )
    },
    logical(1L)
  ))
}

.crs_progress_resolve_message_muffling <- function(state) {
  if (!identical(state$renderer, "single_line")) {
    return(state)
  }

  if (isTRUE(state$message_muffled_checked)) {
    return(state)
  }

  state$message_muffled <- isTRUE(.crs_progress_is_message_muffled())
  state$message_muffled_checked <- TRUE
  state
}

.crs_progress_fmt_num <- function(x) {
  formatC(x, digits = 1L, format = "f")
}

.crs_progress_make_registry <- function() {
  env <- new.env(parent = emptyenv())
  env$next_id <- 0L
  env$active_id <- NULL
  env
}

.crs_progress_registry <- .crs_progress_make_registry()

.crs_progress_reset_registry <- function() {
  .crs_progress_registry$next_id <- 0L
  .crs_progress_registry$active_id <- NULL
  invisible(NULL)
}

.crs_progress_next_id <- function() {
  .crs_progress_registry$next_id <- as.integer(.crs_progress_registry$next_id) + 1L
  sprintf("progress-%d", .crs_progress_registry$next_id)
}

.crs_progress_claim_owner <- function(session_id) {
  active_id <- .crs_progress_registry$active_id
  if (is.null(active_id) || identical(active_id, session_id)) {
    .crs_progress_registry$active_id <- session_id
    return(TRUE)
  }

  FALSE
}

.crs_progress_release_owner <- function(session_id) {
  if (identical(.crs_progress_registry$active_id, session_id)) {
    .crs_progress_registry$active_id <- NULL
  }

  invisible(NULL)
}

.crs_progress_is_rstudio_console <- function() {
  identical(.Platform$GUI, "RStudio") || identical(Sys.getenv("RSTUDIO"), "1")
}

.crs_progress_capability <- function(domain = "general") {
  interactive_console <- isTRUE(.crs_progress_is_interactive())
  rstudio_console <- isTRUE(.crs_progress_is_rstudio_console())

  list(
    domain = domain,
    interactive = interactive_console,
    rstudio = rstudio_console,
    single_line_viable = interactive_console
  )
}

.crs_progress_single_line_surfaces <- function() {
  c("console", "nomad", "bootstrap", "plot", "solver", "cv")
}

.crs_progress_renderer_for_surface <- function(surface, capability) {
  if (surface %in% .crs_progress_single_line_surfaces() &&
      isTRUE(capability$single_line_viable)) {
    return("single_line")
  }

  "legacy"
}

.crs_progress_emit <- function(line) {
  .crs_message(line)
  invisible(NULL)
}

.crs_progress_make_snapshot <- function(state, line, event, now, done = NULL, detail = NULL, render_line = line) {
  list(
    id = state$id,
    surface = state$surface,
    renderer = state$renderer,
    label = state$label,
    kind = if (isTRUE(state$known_total)) "known" else "unknown",
    current = done,
    total = state$total,
    detail = detail,
    line = line,
    render_line = render_line,
    event = event,
    started_at = state$started,
    now = now,
    last_width = state$last_render_width
  )
}

.crs_progress_render_legacy <- function(snapshot, event = c("render", "finish", "abort")) {
  event <- match.arg(event)
  .crs_progress_emit(snapshot$line)
  invisible(snapshot)
}

.crs_progress_single_line_connection <- function() {
  stderr()
}

.crs_progress_has_tty <- function() {
  any(vapply(
    list(stdin(), stdout(), stderr()),
    function(con) {
      tryCatch(isTRUE(isatty(con)), error = function(...) FALSE)
    },
    logical(1L)
  ))
}

.crs_progress_parse_width <- function(value) {
  value <- suppressWarnings(as.integer(value)[1L])
  if (!is.finite(value) || is.na(value) || value < 20L) {
    return(NA_integer_)
  }

  value
}

.crs_progress_terminal_width_probe <- function() {
  if (!isTRUE(.crs_progress_is_interactive()) || isTRUE(.crs_progress_is_rstudio_console())) {
    return(NA_integer_)
  }
  if (!isTRUE(.crs_progress_has_tty())) {
    return(NA_integer_)
  }

  if (identical(.Platform$OS.type, "windows")) {
    return(.crs_progress_parse_width(Sys.getenv("COLUMNS", "")))
  }

  probe <- tryCatch(
    suppressWarnings(system2(
      "sh",
      c("-c", "stty size < /dev/tty 2>/dev/null"),
      stdout = TRUE,
      stderr = TRUE
    )),
    error = function(...) character()
  )
  if (length(probe)) {
    fields <- strsplit(trimws(probe[[1L]]), "[[:space:]]+")[[1L]]
    width <- .crs_progress_parse_width(fields[[length(fields)]])
    if (!is.na(width)) {
      return(width)
    }
  }

  probe <- tryCatch(
    suppressWarnings(system2("tput", "cols", stdout = TRUE, stderr = TRUE)),
    error = function(...) character()
  )
  if (length(probe)) {
    width <- .crs_progress_parse_width(probe[[1L]])
    if (!is.na(width)) {
      return(width)
    }
  }

  .crs_progress_parse_width(Sys.getenv("COLUMNS", ""))
}

.crs_progress_output_width <- function() {
  width <- if (isTRUE(.crs_progress_is_rstudio_console())) {
    .crs_progress_parse_width(getOption("width", 80))
  } else {
    .crs_progress_terminal_width_probe()
  }
  if (is.na(width)) {
    width <- .crs_progress_parse_width(getOption("width", 80))
  }
  if (is.na(width)) {
    width <- 80L
  }

  reserve <- if (isTRUE(.crs_progress_is_rstudio_console())) 4L else 0L
  max(20L, width - reserve)
}

.crs_progress_ellipsize_middle <- function(text, max_width) {
  max_width <- suppressWarnings(as.integer(max_width)[1L])
  if (is.na(max_width) || max_width < 1L) {
    return("")
  }

  if (nchar(text, type = "width") <= max_width) {
    return(text)
  }

  if (max_width <= 3L) {
    return(substr("...", 1L, max_width))
  }

  keep_right <- max(1L, floor((max_width - 3L) / 2L))
  keep_left <- max_width - 3L - keep_right
  paste0(
    substr(text, 1L, keep_left),
    "...",
    substr(text, nchar(text, type = "chars") - keep_right + 1L, nchar(text, type = "chars"))
  )
}

.crs_progress_fit_single_line <- function(line, max_width = .crs_progress_output_width()) {
  if (!is.character(line) || length(line) != 1L || is.na(line)) {
    return(line)
  }

  max_width <- suppressWarnings(as.integer(max_width)[1L])
  if (is.na(max_width) || max_width < 1L) {
    return(line)
  }

  if (nchar(line, type = "width") <= max_width) {
    return(line)
  }

  detail_pos <- regexpr(": ", line, fixed = TRUE)[1L]
  if (detail_pos > 0L) {
    without_detail <- substr(line, 1L, detail_pos - 1L)
    if (nchar(without_detail, type = "width") <= max_width) {
      return(without_detail)
    }
    line <- without_detail
  }

  .crs_progress_ellipsize_middle(line, max_width = max_width)
}

.crs_progress_render_single_line <- function(snapshot, event = c("render", "finish", "abort")) {
  event <- match.arg(event)
  render_line <- snapshot$render_line
  con <- .crs_progress_single_line_connection()
  width <- nchar(render_line, type = "width")

  if (identical(event, "finish")) {
    clear_width <- max(snapshot$last_width, width, .crs_progress_output_width())
    clear_line <- if (clear_width > 0L) strrep(" ", clear_width) else ""
    base::cat("\r", clear_line, "\r\n", file = con, sep = "")
    flush(con)
    flush.console()
    return(invisible(snapshot))
  }

  pad <- max(0L, snapshot$last_width - width)
  suffix <- if (pad > 0L) paste(rep(" ", pad), collapse = "") else ""

  base::cat("\r", render_line, suffix, file = con, sep = "")
  if (identical(event, "abort")) {
    base::cat("\n", file = con, sep = "")
  }
  flush(con)
  flush.console()
  invisible(snapshot)
}

.crs_progress_render <- function(state, line, event, now, done = NULL, detail = NULL) {
  if (!isTRUE(state$enabled) || !isTRUE(state$visible)) {
    return(state)
  }

  state <- .crs_progress_resolve_message_muffling(state)
  if (isTRUE(state$message_muffled)) {
    return(state)
  }

  render_line <- if (identical(state$renderer, "single_line")) {
    .crs_progress_fit_single_line(line)
  } else {
    line
  }

  snapshot <- .crs_progress_make_snapshot(
    state = state,
    line = line,
    render_line = render_line,
    event = event,
    now = now,
    done = done,
    detail = detail
  )

  if (identical(state$renderer, "single_line")) {
    .crs_progress_render_single_line(
      snapshot = snapshot,
      event = if (identical(event, "finish")) "finish" else if (identical(event, "abort")) "abort" else "render"
    )
  } else {
    .crs_progress_render_legacy(
      snapshot = snapshot,
      event = if (identical(event, "finish")) "finish" else if (identical(event, "abort")) "abort" else "render"
    )
  }

  state$rendered <- TRUE
  state$last_render_width <- nchar(render_line, type = "width")
  state$last_line <- line
  state
}

.crs_progress_start_grace_sec <- function(known_total = FALSE, domain = "general") {
  default <- if (identical(domain, "plot")) {
    0.75
  } else if (isTRUE(known_total)) {
    0.75
  } else {
    1.0
  }

  option_name <- if (identical(domain, "plot")) {
    "crs.plot.progress.start.grace.sec"
  } else if (isTRUE(known_total)) {
    "crs.progress.start.grace.known.sec"
  } else {
    "crs.progress.start.grace.unknown.sec"
  }

  val <- suppressWarnings(as.numeric(getOption(option_name, default))[1L])
  if (!is.finite(val) || is.na(val) || val < 0) {
    val <- default
  }

  val
}

.crs_progress_interval_sec <- function(known_total = FALSE, domain = "general") {
  default <- if (identical(domain, "plot")) {
    2.0
  } else if (isTRUE(known_total)) {
    0.5
  } else {
    2.0
  }

  option_name <- if (identical(domain, "plot")) {
    "crs.plot.progress.interval.sec"
  } else if (isTRUE(known_total)) {
    "crs.progress.interval.known.sec"
  } else {
    "crs.progress.interval.unknown.sec"
  }

  val <- suppressWarnings(as.numeric(getOption(option_name, default))[1L])
  if (!is.finite(val) || is.na(val) || val < 0) {
    val <- default
  }

  val
}

.crs_progress_format_known_total <- function(state, done, detail = NULL, now = .crs_progress_now()) {
  total <- state$total
  pct <- if (isTRUE(total > 0)) 100 * done / total else 0
  elapsed <- max(0, now - state$started)
  eta <- if (isTRUE(done > 0) && isTRUE(total >= done)) elapsed * (total - done) / done else 0

  line <- sprintf(
    "%s %s %s/%s (%s%%, elapsed %ss, eta %ss)",
    state$pkg_prefix,
    state$label,
    format(done),
    format(total),
    .crs_progress_fmt_num(pct),
    .crs_progress_fmt_num(elapsed),
    .crs_progress_fmt_num(eta)
  )

  if (!is.null(detail)) {
    line <- paste0(line, ": ", detail)
  }

  line
}

.crs_progress_format_unknown_total <- function(state, done = NULL, detail = NULL, now = .crs_progress_now()) {
  elapsed <- max(0, now - state$started)
  line <- sprintf("%s %s... ", state$pkg_prefix, state$label)

  if (!is.null(done)) {
    line <- paste0(line, "iteration ", format(done), ", ")
  }

  line <- paste0(line, "elapsed ", .crs_progress_fmt_num(elapsed), "s")

  if (!is.null(detail)) {
    line <- paste0(line, ": ", detail)
  }

  line
}

.crs_progress_format_line <- function(state, done = NULL, detail = NULL, now = .crs_progress_now()) {
  if (isTRUE(state$known_total)) {
    .crs_progress_format_known_total(
      state = state,
      done = if (is.null(done)) state$last_done else done,
      detail = detail,
      now = now
    )
  } else {
    .crs_progress_format_unknown_total(
      state = state,
      done = if (is.null(done)) state$last_done else done,
      detail = detail,
      now = now
    )
  }
}

.crs_progress_begin <- function(label, total = NULL, domain = "general", surface = domain) {
  known_total <- !is.null(total)
  throttle_sec <- .crs_progress_interval_sec(known_total = known_total, domain = domain)
  started <- .crs_progress_now()
  capability <- .crs_progress_capability(domain = domain)
  session_id <- .crs_progress_next_id()
  enabled <- .crs_progress_enabled(domain = domain)
  renderer <- .crs_progress_renderer_for_surface(surface = surface, capability = capability)
  visible <- if (!isTRUE(enabled)) {
    FALSE
  } else {
    isTRUE(.crs_progress_claim_owner(session_id))
  }

  list(
    id = session_id,
    enabled = enabled,
    visible = visible,
    pkg_prefix = .crs_io_pkg_prefix(),
    label = label,
    surface = surface,
    total = total,
    known_total = known_total,
    started = started,
    last_emit = started - throttle_sec,
    throttle_sec = throttle_sec,
    last_done = if (known_total) 0 else NULL,
    domain = domain,
    renderer = renderer,
    capability = capability,
    rendered = FALSE,
    start_note = sprintf("%s %s...", .crs_io_pkg_prefix(), label),
    start_note_pending = TRUE,
    start_note_grace_sec = .crs_progress_start_grace_sec(known_total = known_total, domain = domain),
    last_line = NULL,
    last_render_width = 0L,
    last_emitted_done = NULL,
    last_emitted_detail = NULL,
    message_muffled = FALSE,
    message_muffled_checked = FALSE
  )
}

.crs_progress_maybe_emit_start_note <- function(state, now = .crs_progress_now()) {
  if (!isTRUE(state$enabled) || !isTRUE(state$visible) || !isTRUE(state$start_note_pending)) {
    return(state)
  }

  if ((now - state$started) < state$start_note_grace_sec) {
    return(state)
  }

  state <- .crs_progress_render(
    state = state,
    line = state$start_note,
    event = "start",
    now = now
  )
  state$start_note_pending <- FALSE
  state
}

.crs_progress_step_at <- function(state, now, done = NULL, detail = NULL, force = FALSE) {
  if (!is.null(done)) {
    state$last_done <- done
  }

  if (!isTRUE(state$enabled)) {
    return(state)
  }

  state <- .crs_progress_maybe_emit_start_note(state = state, now = now)
  if (isTRUE(state$start_note_pending)) {
    return(state)
  }

  if (!isTRUE(force) && (now - state$last_emit) < state$throttle_sec) {
    return(state)
  }

  line <- .crs_progress_format_line(
    state = state,
    done = if (is.null(done)) state$last_done else done,
    detail = detail,
    now = now
  )

  if (!identical(line, state$last_line)) {
    state <- .crs_progress_render(
      state = state,
      line = line,
      event = "update",
      now = now,
      done = if (is.null(done)) state$last_done else done,
      detail = detail
    )
    state$last_emit <- now
    state$last_emitted_done <- if (is.null(done)) state$last_done else done
    state$last_emitted_detail <- detail
  }

  state
}

.crs_progress_step <- function(state, done = NULL, detail = NULL) {
  .crs_progress_step_at(
    state = state,
    now = .crs_progress_now(),
    done = done,
    detail = detail
  )
}

.crs_progress_end <- function(state, detail = NULL) {
  if (!isTRUE(state$enabled)) {
    .crs_progress_release_owner(state$id)
    return(invisible(state))
  }

  now <- .crs_progress_now()
  if (is.null(state$last_line) &&
      (now - state$started) < state$start_note_grace_sec) {
    .crs_progress_release_owner(state$id)
    return(invisible(state))
  }

  must_clear <- identical(state$renderer, "single_line") && isTRUE(state$rendered)

  if (isTRUE(state$known_total)) {
    done <- if (is.null(state$total)) state$last_done else state$total
    line <- .crs_progress_format_line(state = state, done = done, detail = detail, now = now)
    if (isTRUE(must_clear) ||
        !(identical(done, state$last_emitted_done) && identical(detail, state$last_emitted_detail))) {
      state <- .crs_progress_render(
        state = state,
        line = line,
        event = "finish",
        now = now,
        done = done,
        detail = detail
      )
    }
    .crs_progress_release_owner(state$id)
    return(invisible(state))
  }

  if (!is.null(state$last_done)) {
    line <- .crs_progress_format_line(state = state, done = state$last_done, detail = detail, now = now)
    if (isTRUE(must_clear) ||
        !(identical(state$last_done, state$last_emitted_done) && identical(detail, state$last_emitted_detail))) {
      state <- .crs_progress_render(
        state = state,
        line = line,
        event = "finish",
        now = now,
        done = state$last_done,
        detail = detail
      )
    }
  }

  if (isTRUE(must_clear) && is.null(state$last_done) && !is.null(state$last_line)) {
    state <- .crs_progress_render(
      state = state,
      line = state$last_line,
      event = "finish",
      now = now,
      detail = detail
    )
  }

  .crs_progress_release_owner(state$id)
  invisible(state)
}

.crs_progress_abort <- function(state, detail = NULL) {
  if (!isTRUE(state$enabled)) {
    .crs_progress_release_owner(state$id)
    return(invisible(state))
  }

  now <- .crs_progress_now()
  if (!isTRUE(state$rendered)) {
    .crs_progress_release_owner(state$id)
    return(invisible(state))
  }

  line <- if (!is.null(detail)) paste0(state$pkg_prefix, " ", detail) else state$last_line
  state <- .crs_progress_render(
    state = state,
    line = line,
    event = "abort",
    now = now,
    done = state$last_done,
    detail = detail
  )
  .crs_progress_release_owner(state$id)
  invisible(state)
}

.crs_progress_note <- function(label) {
  if (.crs_progress_enabled()) {
    .crs_progress_emit(sprintf("%s %s", .crs_io_pkg_prefix(), label))
  }

  invisible(NULL)
}

.crs_progress_normalize_manual_line <- function(line) {
  line <- paste(line, collapse = "")
  line <- gsub("\r", "", line, fixed = TRUE)
  line <- gsub("[[:space:]]*\n+[[:space:]]*", " ", line)
  trimws(line, which = "both")
}

.crs_progress_manual_begin <- function(surface = "console") {
  state <- .crs_progress_begin("progress", surface = surface)
  state$start_note_pending <- FALSE
  state$start_note_grace_sec <- 0
  state
}

.crs_progress_status_begin <- function(enabled = TRUE, surface = "console") {
  env <- new.env(parent = emptyenv())
  env$enabled <- isTRUE(enabled)
  env$surface <- surface
  env$state <- NULL
  env
}

.crs_progress_status_update <- function(status, line) {
  if (is.null(status) || !isTRUE(status$enabled)) {
    return(invisible(status))
  }

  if (is.null(status$state)) {
    status$state <- .crs_progress_manual_begin(surface = status$surface)
  }
  status$state <- .crs_progress_manual_update(status$state, line)
  invisible(status)
}

.crs_progress_status_clear <- function(status) {
  if (is.null(status) || !isTRUE(status$enabled)) {
    return(invisible(status))
  }

  if (!is.null(status$state)) {
    .crs_progress_manual_end(status$state)
    status$state <- NULL
  }

  invisible(status)
}

.crs_progress_manual_update <- function(state, line) {
  line <- .crs_progress_normalize_manual_line(line)
  if (!nzchar(line)) {
    return(state)
  }

  if (!isTRUE(state$enabled) || !isTRUE(state$visible)) {
    state$last_line <- .crs_io_prefix_text(line)
    return(state)
  }

  now <- .crs_progress_now()
  line <- .crs_io_prefix_text(line)
  state <- .crs_progress_render(
    state = state,
    line = line,
    event = "update",
    now = now
  )
  state$last_emit <- now
  state$last_emitted_detail <- NULL
  state$last_emitted_done <- NULL
  state$start_note_pending <- FALSE
  state
}

.crs_progress_manual_end <- function(state) {
  if (is.null(state)) {
    return(invisible(NULL))
  }

  if (!isTRUE(state$enabled)) {
    .crs_progress_release_owner(state$id)
    return(invisible(NULL))
  }

  if (isTRUE(state$rendered) && !is.null(state$last_line)) {
    state <- .crs_progress_render(
      state = state,
      line = state$last_line,
      event = "finish",
      now = .crs_progress_now()
    )
  }

  .crs_progress_release_owner(state$id)
  invisible(NULL)
}

.crs_progress_console_runtime <- local({
  env <- new.env(parent = emptyenv())
  env$state <- NULL
  env
})

.crs_progress_console_update <- function(line) {
  state <- .crs_progress_console_runtime$state
  if (is.null(state)) {
    state <- .crs_progress_manual_begin(surface = "console")
  }
  state <- .crs_progress_manual_update(state, line)
  .crs_progress_console_runtime$state <- state
  invisible(NULL)
}

.crs_progress_console_clear <- function() {
  state <- .crs_progress_console_runtime$state
  if (!is.null(state)) {
    .crs_progress_manual_end(state)
    .crs_progress_console_runtime$state <- NULL
  }

  invisible(NULL)
}
