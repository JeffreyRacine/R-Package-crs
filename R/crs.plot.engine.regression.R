## Internal CRS regression plot engines.
##
## The functions here are shadow-only: they are not wired into plot.crs() yet.

.crs_plot_payload_to_legacy_slices <- function(payload) {
  if (!inherits(payload, "crs_plot_payload"))
    stop("payload must inherit from class 'crs_plot_payload'")
  if (!identical(payload$route, "crs") || isTRUE(payload$perspective))
    stop("only non-perspective CRS payloads can be converted to slice data")

  out <- Map(function(slice, nm) {
    if (!(nm %in% names(slice)))
      stop("payload slice does not contain its varying predictor column")
    if (isTRUE(payload$ci) && all(c("lwr", "upr") %in% names(slice))) {
      ans <- data.frame(slice[[nm]], slice$fit, slice$lwr, slice$upr,
                        check.names = FALSE)
      names(ans) <- c(nm, "mean", "lwr", "upr")
      row.names(ans) <- as.character(seq_len(NROW(ans)))
    } else {
      ans <- data.frame(slice[[nm]], slice$fit, check.names = FALSE)
      names(ans) <- c(nm, "mean")
    }
    ans
  }, payload$slices, names(payload$slices))
  names(out) <- names(payload$slices)
  out
}

.crs_plot_slice_ylim <- function(slices, ci = FALSE, common.scale = TRUE) {
  if (!isTRUE(common.scale)) return(NULL)
  vals <- unlist(lapply(slices, function(x) {
    if (isTRUE(ci) && all(c("lwr", "upr") %in% names(x))) {
      as.numeric(unlist(x[, setdiff(names(x), names(x)[1L]), drop = FALSE]))
    } else {
      as.numeric(x[[2L]])
    }
  }), use.names = FALSE)
  range(vals, finite = TRUE)
}

.crs_plot_render_regression_1d <- function(object,
                                           slices,
                                           ci = FALSE,
                                           common.scale = TRUE,
                                           data_overlay = TRUE,
                                           ...) {
  ylim <- .crs_plot_slice_ylim(slices, ci = ci, common.scale = common.scale)
  if (!is.null(object$num.z) || object$num.x > 1L)
    graphics::par(mfrow = grDevices::n2mfrow(length(slices)))

  for (nm in names(slices)) {
    slice <- slices[[nm]]
    x <- slice[[1L]]
    y <- slice[[2L]]
    ylab <- if (is.null(object$tau)) {
      "Conditional Mean"
    } else {
      paste("Conditional Quantile (tau = ", format(object$tau), ")", sep = "")
    }
    local.ylim <- if (is.null(ylim)) {
      range(as.numeric(unlist(slice[, -1L, drop = FALSE])), finite = TRUE)
    } else {
      ylim
    }
    plot.args <- .crs_plot_merge_user_args(
      list(x = x, y = y, xlab = nm, ylab = ylab, ylim = local.ylim,
           type = "l", col = .crs_plot_color("fit"), lwd = 2),
      .crs_plot_user_args(list(...), "plot")
    )
    do.call(graphics::plot, plot.args)

    if (isTRUE(data_overlay) && nm %in% names(object$xz)) {
      yy <- object$y
      xx <- object$xz[[nm]]
      if (!is.factor(xx)) {
        graphics::points(xx, yy, col = .crs_plot_color("data_overlay", 0.25),
                         pch = 16, cex = 0.6)
      }
    }

    if (isTRUE(ci) && all(c("lwr", "upr") %in% names(slice))) {
      graphics::lines(x, slice$lwr, col = .crs_plot_color("interval"),
                      lty = 2)
      graphics::lines(x, slice$upr, col = .crs_plot_color("interval"),
                      lty = 2)
    }
  }

  invisible(slices)
}

.crs_plot_regression_1d_shadow <- function(object,
                                           ...,
                                           .plot_dots_call = NULL) {
  if (is.null(.plot_dots_call))
    .plot_dots_call <- match.call(expand.dots = FALSE)$...
  .crs_plot_validate_public_dots(.plot_dots_call, context = "plot.crs")
  dots <- .crs_plot_normalize_public_dots(list(...), context = "plot.crs")

  plot.behavior <- if (!is.null(dots$plot.behavior)) {
    match.arg(dots$plot.behavior, c("plot", "plot-data", "data"))
  } else {
    "plot"
  }
  ci <- isTRUE(.crs_plot_scalar_default(dots$ci, FALSE))
  num.eval <- as.integer(.crs_plot_scalar_default(dots$num.eval, 100L))
  xtrim <- .crs_plot_scalar_default(dots$xtrim, 0)
  xq <- .crs_plot_scalar_default(dots$xq, 0.5)
  common.scale <- isTRUE(.crs_plot_scalar_default(dots$common.scale, TRUE))
  data_overlay <- isTRUE(.crs_plot_scalar_default(dots$plot.data.overlay, TRUE))
  perspective <- isTRUE(.crs_plot_scalar_default(dots$perspective, FALSE))
  if (isTRUE(perspective))
    stop("modern 2D regression plot route is not implemented yet",
         call. = FALSE)

  payload <- .crs_plot_payload_regression(object = object,
                                          deriv = 0L,
                                          ci = ci,
                                          num.eval = num.eval,
                                          xtrim = xtrim,
                                          xq = xq,
                                          perspective = FALSE,
                                          legacy = FALSE,
                                          display.nomad.progress = FALSE)
  slices <- .crs_plot_payload_to_legacy_slices(payload)

  if (!identical(plot.behavior, "data")) {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar), add = TRUE)
    render.dots <- dots[setdiff(names(dots),
                                c("plot.behavior", "plot.data.overlay",
                                  "num.eval", "xtrim", "xq", "ci",
                                  "common.scale"))]
    do.call(.crs_plot_render_regression_1d,
            c(list(object = object,
                   slices = slices,
                   ci = ci,
                   common.scale = common.scale,
                   data_overlay = data_overlay),
              render.dots))
  }

  if (!identical(plot.behavior, "plot")) return(slices)
  invisible(slices)
}

.crs_plot_regression_1d_public <- function(object,
                                           plot.call,
                                           mean,
                                           deriv,
                                           ci,
                                           plot.errors.method,
                                           num.eval,
                                           xtrim,
                                           xq,
                                           plot.behavior,
                                           common.scale,
                                           persp.rgl,
                                           display.nomad.progress,
                                           display.warnings,
                                           ...) {
  if (!inherits(object, "crs")) stop("object must inherit from class 'crs'")
  if (isTRUE(persp.rgl))
    stop("modern 2D regression plot route is not implemented yet",
         call. = FALSE)
  if (!isTRUE(mean) && identical(deriv, 0))
    mean <- TRUE
  if (!isTRUE(mean) || !identical(as.numeric(deriv), 0))
    stop("plot.view=\"fit\" currently supports fitted mean/quantile plots only",
         call. = FALSE)
  if (isTRUE(ci) && !identical(plot.errors.method, "asymptotic"))
    stop("plot.view=\"fit\" currently supports asymptotic intervals only",
         call. = FALSE)

  dots <- list(...)
  supplied <- names(plot.call)
  raw.dots <- plot.call$...
  dot.names <- names(raw.dots)
  if (is.null(dot.names)) dot.names <- character()
  effective.errors.method <- plot.errors.method
  if ("errors" %in% dot.names) {
    effective.errors.method <- .crs_plot_scalar_match(
      raw.dots$errors,
      c("none", "bootstrap", "asymptotic"),
      "errors"
    )
  }

  bridge <- dots
  if (!("output" %in% dot.names)) bridge$plot.behavior <- plot.behavior
  bridge$ci <- ci
  if (!("neval" %in% dot.names)) bridge$num.eval <- num.eval
  bridge$xtrim <- xtrim
  bridge$xq <- xq
  bridge$common.scale <- common.scale
  bridge$display.nomad.progress <- display.nomad.progress
  bridge$display.warnings <- display.warnings

  if ("plot.behavior" %in% supplied && "output" %in% dot.names)
    bridge$plot.behavior <- plot.behavior
  if (isTRUE(ci) && !identical(effective.errors.method, "asymptotic"))
    stop("plot.view=\"fit\" currently supports asymptotic intervals only",
         call. = FALSE)

  do.call(.crs_plot_regression_1d_shadow,
          c(list(object = object, .plot_dots_call = raw.dots), bridge))
}
