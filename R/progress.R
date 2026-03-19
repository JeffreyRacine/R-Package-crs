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

.crs_progress_iv_enhanced_state <- function(state) {
  isTRUE(state$iv_progress_common)
}

.crs_progress_iv_title <- function() {
  "IV regression"
}

.crs_progress_iv_initialize_state <- function(state) {
  state$iv_progress_common <- TRUE
  state$iv_object_label <- NULL
  state$iv_iteration <- NULL
  state$start_note <- sprintf(
    "%s %s...",
    state$pkg_prefix,
    .crs_progress_iv_title()
  )
  state
}

.crs_progress_iv_format_line <- function(state, now = .crs_progress_now()) {
  elapsed <- max(0, now - state$started)
  fields <- character()

  if (!is.null(state$iv_object_label) && nzchar(state$iv_object_label)) {
    fields <- c(fields, state$iv_object_label)
  }

  if (!is.null(state$iv_iteration)) {
    fields <- c(fields, sprintf("iteration %s", format(state$iv_iteration)))
  }

  fields <- c(fields, sprintf("elapsed %ss", .crs_progress_fmt_num(elapsed)))

  sprintf(
    "%s %s (%s)",
    state$pkg_prefix,
    .crs_progress_iv_title(),
    paste(fields, collapse = ", ")
  )
}

.crs_progress_show_now <- function(state,
                                   done = state$last_done,
                                   detail = state$last_emitted_detail) {
  if (is.null(state)) {
    return(NULL)
  }

  state$last_emit <- -Inf
  .crs_progress_step_at(
    state = state,
    now = .crs_progress_now(),
    done = done,
    detail = detail,
    force = TRUE
  )
}

.crs_plot_progress_enabled <- function() {
  isTRUE(getOption("crs.messages", TRUE)) &&
    isTRUE(getOption("crs.plot.progress", TRUE)) &&
    isTRUE(.crs_progress_is_interactive())
}

.crs_plot_progress_interval_sec <- function() {
  .crs_progress_interval_sec(known_total = FALSE, domain = "plot")
}

.crs_plot_progress_start_grace_sec <- function() {
  .crs_progress_start_grace_sec(known_total = TRUE, domain = "plot")
}

.crs_plot_progress_max_intermediate <- function() {
  val <- suppressWarnings(as.integer(getOption("crs.plot.progress.max.intermediate", 3L))[1L])
  if (is.na(val) || val < 0L) {
    val <- 3L
  }
  val
}

.crs_plot_progress_chunk_cap <- function(total) {
  total <- as.integer(total)
  if (is.na(total) || total < 1L) {
    return(1L)
  }

  max_intermediate <- .crs_plot_progress_max_intermediate()
  if (is.na(max_intermediate) || max_intermediate < 1L) {
    return(total)
  }

  max(1L, as.integer(ceiling(total / (max_intermediate + 1L))))
}

.crs_plot_progress_warmup_max_reps <- function() {
  val <- suppressWarnings(as.integer(getOption("crs.plot.progress.warmup.max.reps", 16L))[1L])
  if (is.na(val) || val < 1L) {
    val <- 16L
  }
  val
}

.crs_plot_progress_warmup_chunk <- function(n,
                                            B,
                                            chunk.size,
                                            progress_enabled = .crs_plot_progress_enabled()) {
  n <- as.integer(n)[1L]
  B <- as.integer(B)[1L]
  chunk.size <- as.integer(chunk.size)[1L]
  if (is.na(n) || n < 1L || is.na(B) || B < 1L || is.na(chunk.size) || chunk.size < 1L) {
    return(1L)
  }
  if (!isTRUE(progress_enabled)) {
    return(min(B, chunk.size))
  }

  warmup.bytes <- 4 * 1024 * 1024
  warmup.chunk <- as.integer(floor(warmup.bytes / (8 * n)))
  if (!is.finite(warmup.chunk) || is.na(warmup.chunk) || warmup.chunk < 1L) {
    warmup.chunk <- 1L
  }
  warmup.chunk <- min(warmup.chunk, .crs_plot_progress_warmup_max_reps())

  min(B, chunk.size, warmup.chunk)
}

.crs_plot_progress_chunk_controller <- function(chunk.size, progress = NULL) {
  chunk.size <- as.integer(chunk.size)[1L]
  if (is.na(chunk.size) || chunk.size < 1L) {
    chunk.size <- 1L
  }

  target.sec <- if (!is.null(progress)) {
    suppressWarnings(as.numeric(progress$throttle_sec)[1L])
  } else {
    NA_real_
  }

  list(
    chunk.size = chunk.size,
    adaptive = isTRUE(.crs_plot_progress_enabled()) &&
      !is.null(progress) &&
      is.finite(target.sec) &&
      !is.na(target.sec) &&
      target.sec > 0,
    target.sec = target.sec
  )
}

.crs_plot_progress_chunk_observe <- function(controller, bsz, elapsed.sec) {
  if (is.null(controller) || !isTRUE(controller$adaptive)) {
    return(controller)
  }

  bsz <- as.integer(bsz)[1L]
  elapsed.sec <- suppressWarnings(as.numeric(elapsed.sec)[1L])
  target.sec <- suppressWarnings(as.numeric(controller$target.sec)[1L])
  if (is.na(bsz) || bsz < 1L || !is.finite(elapsed.sec) || is.na(elapsed.sec) || elapsed.sec <= 0) {
    return(controller)
  }
  if (!is.finite(target.sec) || is.na(target.sec) || target.sec <= 0) {
    return(controller)
  }

  suggested <- as.integer(round(bsz * target.sec / elapsed.sec))
  lower <- max(1L, as.integer(floor(bsz / 4L)))
  upper <- max(lower, as.integer(ceiling(bsz * 4L)))
  if (is.na(suggested) || suggested < 1L) {
    suggested <- lower
  }

  controller$chunk.size <- min(upper, max(lower, suggested))
  controller
}

.crs_plot_progress_checkpoints <- function(total) {
  total <- as.integer(total)
  max_intermediate <- .crs_plot_progress_max_intermediate()
  if (is.na(total) || total < 2L || max_intermediate < 1L) {
    return(integer())
  }

  checkpoints <- unique(as.integer(ceiling(total * seq_len(max_intermediate) / (max_intermediate + 1L))))
  checkpoints[checkpoints >= 1L & checkpoints < total]
}

.crs_plot_progress_begin <- function(total, label) {
  total <- as.integer(total)
  if (is.na(total) || total < 1L || !.crs_plot_progress_enabled()) {
    return(NULL)
  }

  label <- as.character(label)[1L]
  state <- .crs_progress_begin(label = label, total = total, domain = "plot", surface = "plot")
  state$enabled <- isTRUE(.crs_plot_progress_enabled())
  state$throttle_sec <- .crs_plot_progress_interval_sec()
  state$last_emit <- state$started - state$throttle_sec
  state$start_note_grace_sec <- .crs_plot_progress_start_grace_sec()
  state$start_note_consumes_throttle <- TRUE
  state$checkpoints <- .crs_plot_progress_checkpoints(total = total)
  state$next_checkpoint_idx <- 1L
  state
}

.crs_plot_stage_progress_begin <- function(total, label) {
  state <- .crs_plot_progress_begin(total = total, label = label)
  if (is.null(state)) {
    return(NULL)
  }

  .crs_progress_show_now(state = state, done = 0L)
}

.crs_plot_progress_tick <- function(state, done, detail = NULL, force = FALSE) {
  if (is.null(state)) {
    return(NULL)
  }

  done <- as.integer(done)
  if (is.na(done)) {
    done <- 0L
  }
  done <- max(0L, min(state$total, done))
  state$last_done <- done

  now <- .crs_progress_now()
  state <- .crs_progress_maybe_emit_start_note(state = state, now = now)
  emit_now <- isTRUE(force)

  if (!isTRUE(force)) {
    checkpoints <- state$checkpoints
    next_idx <- state$next_checkpoint_idx
    reached_checkpoint <- FALSE
    if (length(checkpoints) && !is.na(next_idx) && next_idx <= length(checkpoints)) {
      next_checkpoint <- checkpoints[[next_idx]]
      if (done >= next_checkpoint) {
        reached <- max(which(checkpoints <= done))
        state$next_checkpoint_idx <- as.integer(reached + 1L)
        reached_checkpoint <- TRUE
      }
    }

    emitted_done <- if (is.null(state$last_emitted_done)) 0L else as.integer(state$last_emitted_done)
    advanced <- isTRUE(done > emitted_done)
    time_ready <- isTRUE(advanced) &&
      !isTRUE(state$start_note_pending) &&
      ((now - state$last_emit) >= state$throttle_sec)

    if (!isTRUE(reached_checkpoint) && !isTRUE(time_ready)) {
      return(state)
    }

    emit_now <- isTRUE(reached_checkpoint) || isTRUE(time_ready)
  }

  .crs_progress_step_at(
    state = state,
    now = now,
    done = done,
    detail = detail,
    force = emit_now
  )
}

.crs_plot_progress_end <- function(state, detail = NULL) {
  if (is.null(state)) {
    return(invisible(NULL))
  }

  .crs_progress_end(state = state, detail = detail)
  invisible(NULL)
}

.crs_plot_bootstrap_stage_label <- function(stage,
                                            method_label = NULL,
                                            target_label = NULL) {
  stage <- as.character(stage)[1L]
  method_label <- if (is.null(method_label)) NULL else as.character(method_label)[1L]
  target_label <- if (is.null(target_label)) NULL else as.character(target_label)[1L]

  base <- if (!is.null(method_label) && nzchar(method_label)) {
    sprintf("%s %s", stage, method_label)
  } else {
    stage
  }

  if (!is.null(target_label) && nzchar(target_label)) {
    sprintf("%s (%s)", base, target_label)
  } else {
    base
  }
}

.crs_plot_progress_target_name <- function(name, fallback) {
  name <- if (is.null(name)) NULL else as.character(name)[1L]
  if (is.null(name) || !nzchar(name) || is.na(name)) {
    fallback
  } else {
    name
  }
}

.crs_plot_progress_target_label <- function(target_name = NULL,
                                            index = 1L,
                                            total = 1L) {
  index <- suppressWarnings(as.integer(index)[1L])
  total <- suppressWarnings(as.integer(total)[1L])
  target_name <- if (is.null(target_name)) NULL else as.character(target_name)[1L]

  if (is.na(total) || total < 1L) {
    total <- 1L
  }
  if (is.na(index) || index < 1L) {
    index <- 1L
  }

  if (total <= 1L) {
    return(NULL)
  }

  if (!is.null(target_name) && nzchar(target_name) && !is.na(target_name)) {
    sprintf("%s %d/%d", target_name, index, total)
  } else {
    sprintf("surf %d/%d", index, total)
  }
}

.crs_plot_regression_bootstrap_target_label <- function(object,
                                                        slice.index,
                                                        gradients = FALSE) {
  slice.index <- suppressWarnings(as.integer(slice.index)[1L])
  total <- suppressWarnings(as.integer(NCOL(object$xz))[1L])
  if (is.na(total) || total < 1L) {
    total <- 1L
  }

  if (is.na(slice.index) || slice.index <= 0L) {
    return(.crs_plot_progress_target_label(index = 1L, total = total))
  }

  target_name <- .crs_plot_progress_target_name(
    if (!is.null(names(object$xz)) && length(names(object$xz)) >= slice.index) names(object$xz)[[slice.index]] else NULL,
    sprintf("x%d", slice.index)
  )
  if (isTRUE(gradients)) {
    target_name <- sprintf("grad %s", target_name)
  }

  .crs_plot_progress_target_label(
    target_name = target_name,
    index = slice.index,
    total = total
  )
}

.crs_plot_glp_bootstrap_target_label <- function(object,
                                                 slice.index,
                                                 gradients = FALSE) {
  slice.index <- suppressWarnings(as.integer(slice.index)[1L])
  total <- suppressWarnings(as.integer(NCOL(object$x))[1L])
  if (is.na(total) || total < 1L) {
    total <- 1L
  }

  if (is.na(slice.index) || slice.index <= 0L) {
    return(.crs_plot_progress_target_label(index = 1L, total = total))
  }

  target_name <- .crs_plot_progress_target_name(
    if (!is.null(object$xnames) && length(object$xnames) >= slice.index) object$xnames[[slice.index]] else NULL,
    sprintf("x%d", slice.index)
  )
  if (isTRUE(gradients)) {
    target_name <- sprintf("grad %s", target_name)
  }

  .crs_plot_progress_target_label(
    target_name = target_name,
    index = slice.index,
    total = total
  )
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
  if (.crs_progress_iv_enhanced_state(state)) {
    return(.crs_progress_iv_format_line(
      state = state,
      now = now
    ))
  }

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
