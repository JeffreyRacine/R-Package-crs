.crs_nomad_scalar <- function(x) {
  if (is.null(x) || !length(x)) return(NA_real_)
  suppressWarnings(as.numeric(x)[1L])
}

.crs_nomad_is_missing_number <- function(x) {
  length(x) != 1L || is.na(x)
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
    label <- paste("\nNumber of Function Evaluations: ", format(nfe), sep = "")
    parts <- character()
    if (!.crs_nomad_is_missing_number(summary$total.evaluations)) {
      parts <- c(parts, paste("total = ", format(summary$total.evaluations), sep = ""))
    }
    if (!.crs_nomad_is_missing_number(summary$cache.hits)) {
      parts <- c(parts, paste("cache hits = ", format(summary$cache.hits), sep = ""))
    }
    if (!.crs_nomad_is_missing_number(summary$callback.evaluations)) {
      parts <- c(parts, paste("callbacks = ", format(summary$callback.evaluations), sep = ""))
    }
    if (length(parts)) {
      label <- paste0(label, " (", paste(parts, collapse = ", "), ")")
    }
    cat(label)
  }

  cache.parts <- character()
  if (!.crs_nomad_is_missing_number(summary$cache.size)) {
    cache.parts <- c(cache.parts, paste("native size = ", format(summary$cache.size), sep = ""))
  }
  if (!.crs_nomad_is_missing_number(summary$callback.cache.hits)) {
    cache.parts <- c(cache.parts, paste("callback hits = ", format(summary$callback.cache.hits), sep = ""))
  }
  if (!.crs_nomad_is_missing_number(summary$callback.cache.misses)) {
    cache.parts <- c(cache.parts, paste("callback misses = ", format(summary$callback.cache.misses), sep = ""))
  }
  if (length(cache.parts)) {
    cat(paste("\nNOMAD Cache Summary: ", paste(cache.parts, collapse = ", "), sep = ""))
  }

  invisible(TRUE)
}
