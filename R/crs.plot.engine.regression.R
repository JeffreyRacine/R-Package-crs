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

.crs_plot_payload_to_legacy_surface <- function(payload) {
  if (!inherits(payload, "crs_plot_payload"))
    stop("payload must inherit from class 'crs_plot_payload'")
  if (!identical(payload$route, "crs") || !isTRUE(payload$perspective))
    stop("only perspective CRS payloads can be converted to surface data")

  xnames <- names(payload$data)[seq_len(2L)]
  out <- data.frame(payload$data[, xnames, drop = FALSE], payload$z)
  list(out)
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

.crs_plot_render_regression_surface <- function(object,
                                                payload,
                                                renderer = c("base", "rgl"),
                                                data_overlay = TRUE,
                                                ...) {
  renderer <- match.arg(renderer)
  zlab <- if (is.null(object$tau)) {
    "Conditional Mean"
  } else {
    paste("Conditional Quantile (tau = ", format(object$tau), ")", sep = "")
  }
  main <- zlab

  if (identical(renderer, "rgl")) {
    col <- .crs_plot_rgl_surface_colors(payload$z)
    return(.crs_plot_render_surface_rgl(
      x = payload$x,
      y = payload$y,
      z = payload$z,
      xlab = names(object$xz)[1L],
      ylab = names(object$xz)[2L],
      zlab = "Y",
      main = main,
      col = col
    ))
  }

  col <- matrix(.crs_plot_surface_colors(payload$z), nrow = NROW(payload$z))
  persp.args <- .crs_plot_merge_user_args(
    list(x = payload$x,
         y = payload$y,
         z = payload$z,
         xlab = names(object$xz)[1L],
         ylab = names(object$xz)[2L],
         zlab = "Y",
         main = main,
         col = col,
         border = .crs_plot_color("surface_border"),
         ticktype = "detailed",
         theta = 45,
         phi = 30),
    .crs_plot_user_args(list(...), "persp")
  )
  persp.mat <- do.call(graphics::persp, persp.args)

  if (isTRUE(data_overlay)) {
    dots <- list(...)
    points.args <- .crs_plot_merge_user_args(
      list(x = object$xz[, 1L],
           y = object$xz[, 2L],
           z = object$y,
           pmat = persp.mat,
           col = .crs_plot_color("data_overlay", 0.35),
           pch = 16,
           cex = 0.5),
      .crs_plot_user_args(dots, "points")
    )
    xy <- do.call(grDevices::trans3d, points.args[c("x", "y", "z", "pmat")])
    points.args <- points.args[setdiff(names(points.args), c("x", "y", "z", "pmat"))]
    do.call(graphics::points, c(list(x = xy$x, y = xy$y), points.args))
  }

  invisible(persp.mat)
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

.crs_plot_regression_surface_shadow <- function(object,
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
  if (isTRUE(ci))
    stop("modern 2D regression plot route does not yet support intervals",
         call. = FALSE)
  num.eval <- as.integer(.crs_plot_scalar_default(dots$num.eval, 100L))
  xtrim <- .crs_plot_scalar_default(dots$xtrim, 0)
  renderer <- .crs_plot_scalar_default(dots$renderer, "base")
  renderer <- match.arg(renderer, c("base", "rgl"))
  data_overlay <- isTRUE(.crs_plot_scalar_default(dots$plot.data.overlay, TRUE))

  payload <- .crs_plot_payload_regression(object = object,
                                          deriv = 0L,
                                          ci = FALSE,
                                          num.eval = num.eval,
                                          xtrim = xtrim,
                                          perspective = TRUE,
                                          legacy = FALSE,
                                          display.nomad.progress = FALSE)
  surface <- .crs_plot_payload_to_legacy_surface(payload)

  if (!identical(plot.behavior, "data")) {
    render.dots <- dots[setdiff(names(dots),
                                c("plot.behavior", "plot.data.overlay",
                                  "num.eval", "xtrim", "ci", "perspective",
                                  "renderer"))]
    do.call(.crs_plot_render_regression_surface,
            c(list(object = object,
                   payload = payload,
                   renderer = renderer,
                   data_overlay = data_overlay),
              render.dots))
  }

  if (!identical(plot.behavior, "plot")) return(surface)
  invisible(surface)
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
  surface.request <- isTRUE(persp.rgl)
  if ("perspective" %in% dot.names)
    surface.request <- isTRUE(raw.dots$perspective)
  if ("renderer" %in% dot.names)
    surface.request <- TRUE
  if (isTRUE(persp.rgl) && "renderer" %in% dot.names &&
      !identical(as.character(raw.dots$renderer), "rgl"))
    stop("cannot supply persp.rgl=TRUE with renderer other than \"rgl\"",
         call. = FALSE)
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
  if (isTRUE(persp.rgl) && !("renderer" %in% dot.names))
    bridge$renderer <- "rgl"

  if (isTRUE(surface.request)) {
    if (isTRUE(ci))
      stop("plot.view=\"fit\" surface route currently supports point estimates only",
           call. = FALSE)
    return(do.call(.crs_plot_regression_surface_shadow,
                   c(list(object = object, .plot_dots_call = raw.dots),
                     bridge)))
  }

  do.call(.crs_plot_regression_1d_shadow,
          c(list(object = object, .plot_dots_call = raw.dots), bridge))
}
