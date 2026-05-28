.crs_nomad_clamp_to_bounds <- function(x, lb, ub) {
  if (is.finite(lb) && x < lb) {
    x <- lb
  }
  if (is.finite(ub) && x > ub) {
    x <- ub
  }
  x
}

.crs_nomad_coerce_start_value <- function(x, type, lb, ub) {
  if (identical(type, 3L)) {
    x <- if (x >= 0.5) 1 else 0
  } else if (type %in% c(1L, 2L)) {
    x <- round(x)
  }

  .crs_nomad_clamp_to_bounds(x, lb, ub)
}

.crs_nomad_build_start_matrix_fallback <- function(x0,
                                                   nstart,
                                                   bbin,
                                                   lb,
                                                   ub) {
  n <- length(bbin)
  starts <- matrix(0, nrow = nstart, ncol = n)

  for (j in seq_len(nstart)) {
    for (i in seq_len(n)) {
      lo <- if (is.finite(lb[i])) lb[i] else -1
      hi <- if (is.finite(ub[i])) ub[i] else 1
      if (hi < lo) {
        tmp <- lo
        lo <- hi
        hi <- tmp
      }
      starts[j, i] <- .crs_nomad_coerce_start_value(
        x = runif(1L, min = lo, max = hi),
        type = as.integer(bbin[i]),
        lb = lb[i],
        ub = ub[i]
      )
    }
  }

  if (is.null(x0)) {
    return(starts)
  }

  x0 <- as.numeric(x0)
  if (length(x0) == n) {
    starts[1L, ] <- vapply(
      seq_len(n),
      function(i) {
        .crs_nomad_coerce_start_value(
          x = x0[i],
          type = as.integer(bbin[i]),
          lb = lb[i],
          ub = ub[i]
        )
      },
      numeric(1L)
    )
    return(starts)
  }

  if (length(x0) >= n * nstart) {
    x0.mat <- matrix(x0[seq_len(n * nstart)], nrow = nstart, ncol = n)
    for (j in seq_len(nstart)) {
      starts[j, ] <- vapply(
        seq_len(n),
        function(i) {
          .crs_nomad_coerce_start_value(
            x = x0.mat[j, i],
            type = as.integer(bbin[i]),
            lb = lb[i],
            ub = ub[i]
          )
        },
        numeric(1L)
      )
    }
  }

  starts
}

.crs_nomad_capture_start_matrix <- function(x0,
                                            nstart,
                                            bbin,
                                            bbout,
                                            lb,
                                            ub,
                                            random.seed,
                                            opts) {
  nstart <- max(1L, as.integer(nstart)[1L])
  n <- length(bbin)

  if (nstart <= 1L) {
    return(matrix(as.numeric(x0), nrow = 1L, ncol = n, byrow = TRUE))
  }

  seed.state <- .crs_capture_seed()
  on.exit(.crs_restore_seed(seed.state), add = TRUE)

  probe.opts <- opts
  probe.opts$DISPLAY_DEGREE <- 0L
  probe.opts$MAX_BB_EVAL <- 1L

  probe.env <- new.env(parent = baseenv())
  probe.env$captured_points <- list()
  probe.env$eval.capture <- function(x) {
      captured_points[[length(captured_points) + 1L]] <- as.numeric(x)
      0
  }

  try(
    suppressWarnings(
      snomadr(
        eval.f = probe.env$eval.capture,
        n = n,
        x0 = as.numeric(x0),
        bbin = bbin,
        bbout = bbout,
        lb = lb,
        ub = ub,
        nmulti = as.integer(nstart),
        random.seed = random.seed,
        opts = probe.opts,
        display.nomad.progress = FALSE,
        snomadr.environment = probe.env
      )
    ),
    silent = TRUE
  )

  if (length(probe.env$captured_points) >= nstart) {
    starts <- matrix(
      unlist(probe.env$captured_points[seq_len(nstart)], use.names = FALSE),
      nrow = nstart,
      ncol = n,
      byrow = TRUE
    )
    return(starts)
  }

  if (!is.null(random.seed) && length(random.seed)) {
    set.seed(as.integer(random.seed[1L]))
  }
  .crs_nomad_build_start_matrix_fallback(
    x0 = x0,
    nstart = nstart,
    bbin = bbin,
    lb = lb,
    ub = ub
  )
}

.crs_nomad_format_int_vec <- function(x) {
  if (is.null(x) || !length(x)) {
    return("pending")
  }

  sprintf("(%s)", paste(format(as.integer(x)), collapse = ","))
}

.crs_nomad_format_num_vec <- function(x) {
  if (is.null(x) || !length(x)) {
    return("pending")
  }

  sprintf("(%s)", paste(vapply(as.numeric(x), format, character(1L)), collapse = ","))
}

.crs_nomad_progress_show_start <- function() {
  isTRUE(getOption("crs.nomad.progress.show.start", TRUE))
}

.crs_nomad_start_equal <- function(x,
                                   y,
                                   tol = sqrt(.Machine$double.eps)) {
  x <- as.numeric(x)
  y <- as.numeric(y)

  isTRUE(all.equal(x, y, tolerance = tol, check.attributes = FALSE))
}

.crs_nomad_progress_begin <- function(progress.status,
                                      label,
                                      starts = NULL,
                                      aux.label = NULL,
                                      aux.type = c("int", "num")) {
  aux.type <- match.arg(aux.type)
  env <- new.env(parent = emptyenv())
  env$enabled <- !is.null(progress.status) && isTRUE(progress.status$enabled)
  env$progress.status <- progress.status
  env$label <- label
  env$starts <- starts
  env$nmulti <- if (is.null(starts)) 1L else nrow(starts)
  env$next.restart <- 1L
  env$restart.index <- 1L
  env$restart.eval <- 0L
  env$iteration <- 0L
  env$started <- .crs_progress_now()
  env$best.objective <- Inf
  env$best.degree <- NULL
  env$best.segments <- NULL
  env$best.aux <- NULL
  env$aux.label <- aux.label
  env$aux.type <- aux.type
  env
}

.crs_nomad_progress_format_aux <- function(progress,
                                           x) {
  if (is.null(progress$aux.label) || is.null(x) || !length(x)) {
    return(NULL)
  }

  formatter <- if (identical(progress$aux.type, "num")) {
    .crs_nomad_format_num_vec
  } else {
    .crs_nomad_format_int_vec
  }

  sprintf("%s %s", progress$aux.label, formatter(x))
}

.crs_nomad_progress_update_line <- function(progress,
                                            degree,
                                            segments,
                                            aux = NULL,
                                            value = NULL,
                                            start = FALSE) {
  if (is.null(progress) || !isTRUE(progress$enabled)) {
    return(invisible(NULL))
  }

  elapsed <- max(0, .crs_progress_now() - progress$started)
  fields <- c(
    sprintf("multistart %s/%s", format(progress$restart.index), format(progress$nmulti)),
    sprintf("eval %s", format(progress$restart.eval)),
    sprintf("deg %s", .crs_nomad_format_int_vec(degree)),
    sprintf("seg %s", .crs_nomad_format_int_vec(segments))
  )

  aux.field <- .crs_nomad_progress_format_aux(progress, aux)
  if (!is.null(aux.field)) {
    fields <- c(fields, aux.field)
  }

  fields <- c(
    fields,
    sprintf("best deg %s", .crs_nomad_format_int_vec(progress$best.degree)),
    sprintf("best seg %s", .crs_nomad_format_int_vec(progress$best.segments))
  )

  best.aux.field <- .crs_nomad_progress_format_aux(progress, progress$best.aux)
  if (!is.null(best.aux.field)) {
    fields <- c(fields, paste("best", best.aux.field))
  }

  if (!is.null(value) && is.finite(value)) {
    fields <- c(fields, sprintf("fv=%s", format(value)))
  } else if (isTRUE(start)) {
    fields <- c(fields, "fv=pending")
  }

  .crs_progress_status_update(
    progress$progress.status,
    sprintf(
      "%s... iteration %s, elapsed %ss: %s",
      progress$label,
      format(progress$iteration),
      .crs_progress_fmt_num(elapsed),
      paste(fields, collapse = ", ")
    )
  )

  invisible(NULL)
}

.crs_nomad_progress_maybe_restart <- function(progress,
                                              input,
                                              degree,
                                              segments,
                                              aux = NULL) {
  if (is.null(progress) || !isTRUE(progress$enabled)) {
    return(invisible(NULL))
  }

  idx <- progress$next.restart
  if (is.null(progress$starts) || idx > nrow(progress$starts)) {
    if (progress$iteration == 0L && progress$restart.eval == 0L) {
      if (isTRUE(.crs_nomad_progress_show_start())) {
        progress$iteration <- 1L
        .crs_nomad_progress_update_line(
          progress = progress,
          degree = degree,
          segments = segments,
          aux = aux,
          value = NULL,
          start = TRUE
        )
      }
    }
    return(invisible(NULL))
  }

  if (.crs_nomad_start_equal(input, progress$starts[idx, ])) {
    progress$restart.index <- idx
    progress$restart.eval <- 0L
    progress$next.restart <- idx + 1L
    if (isTRUE(.crs_nomad_progress_show_start())) {
      progress$iteration <- progress$iteration + 1L
      .crs_nomad_progress_update_line(
        progress = progress,
        degree = degree,
        segments = segments,
        aux = aux,
        value = NULL,
        start = TRUE
      )
    }
  }

  invisible(NULL)
}

.crs_nomad_progress_eval <- function(progress,
                                     degree,
                                     segments,
                                     aux = NULL,
                                     value) {
  if (is.null(progress) || !isTRUE(progress$enabled)) {
    return(invisible(NULL))
  }

  progress$iteration <- progress$iteration + 1L
  progress$restart.eval <- progress$restart.eval + 1L

  if (!is.null(value) && is.finite(value) &&
      (is.null(progress$best.degree) || isTRUE(value < progress$best.objective))) {
    progress$best.objective <- value
    progress$best.degree <- as.integer(degree)
    progress$best.segments <- as.integer(segments)
    progress$best.aux <- aux
  }

  .crs_nomad_progress_update_line(
    progress = progress,
    degree = degree,
    segments = segments,
    aux = aux,
    value = value
  )

  invisible(NULL)
}

.crs_nomad_progress_finish <- function(progress) {
  if (is.null(progress) || !isTRUE(progress$enabled)) {
    return(invisible(NULL))
  }

  elapsed <- max(0, .crs_progress_now() - progress$started)
  fields <- c(
    sprintf("best deg %s", .crs_nomad_format_int_vec(progress$best.degree)),
    sprintf("best seg %s", .crs_nomad_format_int_vec(progress$best.segments))
  )

  best.aux.field <- .crs_nomad_progress_format_aux(progress, progress$best.aux)
  if (!is.null(best.aux.field)) {
    fields <- c(fields, paste("best", best.aux.field))
  }

  fields <- c(
    fields,
    sprintf(
      "fv=%s",
      if (is.finite(progress$best.objective)) format(progress$best.objective) else "pending"
    )
  )

  .crs_progress_status_update(
    progress$progress.status,
    sprintf(
      "%s... elapsed %ss: %s",
      progress$label,
      .crs_progress_fmt_num(elapsed),
      paste(fields, collapse = ", ")
    )
  )

  progress$enabled <- FALSE
  invisible(NULL)
}
