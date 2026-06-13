.crs_nomad_scalar <- function(x) {
  if (is.null(x) || !length(x)) return(NA_real_)
  suppressWarnings(as.numeric(x)[1L])
}

.crs_nomad_is_missing_number <- function(x) {
  length(x) != 1L || is.na(x)
}

.crs_nomad_format_count <- function(x) {
  format(round(x), big.mark = ",", trim = TRUE, scientific = FALSE)
}

.crs_nomad_format_percent <- function(hits, requests) {
  format(round(100 * as.numeric(hits)[1L] / as.numeric(requests)[1L], 1L),
         nsmall = 1L, trim = TRUE)
}

.crs_nomad_native_cache_line <- function(summary) {
  hits <- summary$cache.hits
  requests <- summary$total.evaluations

  if (!.crs_nomad_is_missing_number(hits) &&
      !.crs_nomad_is_missing_number(requests) &&
      is.finite(hits) &&
      is.finite(requests) &&
      requests > 0 &&
      hits > 0 &&
      requests >= hits) {
    return(paste0(
      "NOMAD cache: ",
      .crs_nomad_format_count(hits),
      " repeated point lookups avoided out of ",
      .crs_nomad_format_count(requests),
      " (",
      .crs_nomad_format_percent(hits, requests),
      "%)"
    ))
  }

  if (!.crs_nomad_is_missing_number(hits) &&
      is.finite(hits) &&
      hits > 0) {
    return(paste0("NOMAD cache: ", .crs_nomad_format_count(hits),
                  " repeated point lookups avoided"))
  }

  NULL
}

.crs_nomad_restart_cache_line <- function(summary) {
  hits <- summary$callback.cache.hits
  callbacks <- summary$callback.evaluations

  if (.crs_nomad_is_missing_number(hits) ||
      !is.finite(hits) ||
      hits <= 0) {
    return(NULL)
  }

  if (!.crs_nomad_is_missing_number(callbacks) &&
      is.finite(callbacks) &&
      callbacks > 0 &&
      callbacks >= hits) {
    return(paste0(
      "CRS restart cache: ",
      .crs_nomad_format_count(hits),
      " repeated objective callbacks avoided out of ",
      .crs_nomad_format_count(callbacks),
      " (",
      .crs_nomad_format_percent(hits, callbacks),
      "%)"
    ))
  }

  paste0("CRS restart cache: ", .crs_nomad_format_count(hits),
         " repeated objective callbacks avoided")
}

.crs_nomad_summary_from_solution <- function(x) {
  if (is.null(x)) return(NULL)

  out <- list(
    best.restart = .crs_nomad_scalar(x$best.restart),
    num.feval = .crs_nomad_scalar(x$bbe),
    cache.hits = .crs_nomad_scalar(x$cache.hits),
    cache.size = .crs_nomad_scalar(x$cache.size),
    callback.evaluations = .crs_nomad_scalar(x$callback.evaluations),
    total.evaluations = .crs_nomad_scalar(x$total.evaluations),
    callback.cache.enabled = if (is.null(x$restart.callback.cache.enabled)) NA
      else isTRUE(x$restart.callback.cache.enabled),
    callback.cache.hits = .crs_nomad_scalar(x$restart.callback.cache.hits),
    callback.cache.misses = .crs_nomad_scalar(x$restart.callback.cache.misses),
    callback.cache.size = .crs_nomad_scalar(x$restart.callback.cache.size)
  )

  out
}

.crs_nomad_summary_from_cv <- function(x) {
  if (!is.null(x) && !is.null(x$nomad.summary)) return(x$nomad.summary)
  if (is.null(x) || is.null(x$nomad.restart.evaluations)) return(NULL)
  evals <- x$nomad.restart.evaluations
  out <- list(
    best.restart = if (is.null(x$nomad.best.restart)) NA_integer_
      else as.integer(x$nomad.best.restart)[1L],
    num.feval = NA_real_,
    cache.hits = NA_real_,
    cache.size = NA_real_,
    callback.evaluations = NA_real_,
    total.evaluations = NA_real_,
    callback.cache.enabled = NA,
    callback.cache.hits = NA_real_,
    callback.cache.misses = NA_real_,
    callback.cache.size = NA_real_
  )

  if (is.data.frame(evals)) {
    sum_col <- function(nm) {
      if (!nm %in% names(evals)) return(NA_real_)
      vals <- suppressWarnings(as.numeric(evals[[nm]]))
      if (!length(vals) || all(is.na(vals))) NA_real_ else sum(vals, na.rm = TRUE)
    }
    last_col <- function(nm) {
      if (!nm %in% names(evals)) return(NA_real_)
      vals <- suppressWarnings(as.numeric(evals[[nm]]))
      vals <- vals[!is.na(vals)]
      if (!length(vals)) NA_real_ else vals[[length(vals)]]
    }
    out$num.feval <- sum_col("bbe")
    out$cache.hits <- sum_col("cache.hits")
    out$callback.evaluations <- sum_col("callback.evaluations")
    out$total.evaluations <- sum_col("total.evaluations")
    out$cache.size <- last_col("cache.size")
    out$callback.cache.hits <- sum_col("callback.cache.hits")
    out$callback.cache.misses <- sum_col("callback.cache.misses")
    out$callback.cache.size <- sum_col("callback.cache.entries.added")
  }

  out
}

.crs_nomad_attach_summary <- function(object, cv.return = NULL) {
  summary <- .crs_nomad_summary_from_cv(cv.return)
  if (!is.null(summary)) {
    object$nomad.summary <- summary
  }
  object
}

.crs_nomad_summary_print <- function(object) {
  summary <- object$nomad.summary
  if (is.null(summary)) return(invisible(FALSE))

  nfe <- summary$num.feval
  if (!.crs_nomad_is_missing_number(nfe)) {
    cat(paste("\nNumber of Function Evaluations: ",
              .crs_nomad_format_count(nfe),
              sep = ""))
  }

  native.cache.line <- .crs_nomad_native_cache_line(summary)
  if (!is.null(native.cache.line)) {
    cat(paste("\n", native.cache.line, sep = ""))
  }

  restart.cache.line <- .crs_nomad_restart_cache_line(summary)
  if (!is.null(restart.cache.line)) {
    cat(paste("\n", restart.cache.line, sep = ""))
  }

  invisible(TRUE)
}
