## Internal CRS curve plot engines for IV and CLSD opt-in routes.

.crs_plot_curve_behavior <- function(dots) {
  if (!is.null(dots$output)) {
    .crs_plot_scalar_match(dots$output, c("plot", "plot-data", "data"),
                           "output")
  } else if (!is.null(dots$plot.behavior)) {
    .crs_plot_scalar_match(dots$plot.behavior, c("plot", "plot-data", "data"),
                           "plot.behavior")
  } else {
    "plot"
  }
}

.crs_plot_curve_draw <- function(frame,
                                 xname,
                                 yname,
                                 ci = FALSE,
                                 data = NULL,
                                 data_overlay = FALSE,
                                 ...) {
  dots <- list(...)
  plot.args <- .crs_plot_merge_user_args(
    list(x = frame[[xname]],
         y = frame$fit,
         type = "l",
         xlab = xname,
         ylab = yname,
         lwd = 2,
         col = .crs_plot_color("fit")),
    .crs_plot_user_args(dots, "plot")
  )
  if (isTRUE(data_overlay) && !is.null(data)) {
    plot.args$ylim <- range(c(data$y, frame$fit, frame$lwr, frame$upr),
                            finite = TRUE)
    plot.args$type <- "n"
  }
  do.call(graphics::plot, plot.args)

  if (isTRUE(data_overlay) && !is.null(data)) {
    points.args <- .crs_plot_merge_user_args(
      list(x = data[[xname]],
           y = data$y,
           col = .crs_plot_color("data_overlay", 0.30),
           pch = 16,
           cex = 0.6),
      .crs_plot_user_args(dots, "points")
    )
    do.call(graphics::points, points.args)
    lines.args <- .crs_plot_merge_user_args(
      list(x = frame[[xname]],
           y = frame$fit,
           col = .crs_plot_color("fit"),
           lwd = 2),
      .crs_plot_user_args(dots, "lines")
    )
    do.call(graphics::lines, lines.args)
  }

  if (isTRUE(ci) && all(c("lwr", "upr") %in% names(frame))) {
    graphics::lines(frame[[xname]], frame$lwr,
                    col = .crs_plot_color("interval"), lty = 2)
    graphics::lines(frame[[xname]], frame$upr,
                    col = .crs_plot_color("interval"), lty = 2)
  }

  invisible(frame)
}

.crs_plot_iv_public <- function(object,
                                plot.call,
                                plot.data = FALSE,
                                ci = FALSE,
                                deriv = FALSE,
                                xtrim = 0,
                                ...) {
  raw.dots <- plot.call$...
  .crs_plot_validate_public_dots(raw.dots, context = "plot.crsiv",
                                 extra = c("plot.data", "line.lwd",
                                           "line.lty", "line.col", "ci.lty",
                                           "ci.lwd", "ci.col"))
  dots <- list(...)
  .crs_plot_reject_curve_controls(dots, context = "plot.crsiv",
                                  allow.data.overlay = TRUE)
  plot.behavior <- .crs_plot_curve_behavior(dots)
  if (!is.null(dots$output)) dots$output <- NULL
  if (!is.null(dots$plot.behavior)) dots$plot.behavior <- NULL
  data_overlay <- if (!is.null(dots$data_overlay)) {
    .crs_plot_match_flag(dots$data_overlay, "data_overlay")
  } else {
    plot.data
  }
  dots$data_overlay <- NULL

  payload <- .crs_plot_payload_iv(object, deriv = deriv, ci = ci,
                                  xtrim = xtrim)
  frame <- payload$data
  names(frame)[1L] <- payload$xname
  yname <- if (isTRUE(deriv)) paste("d", "y", "/d", payload$xname, sep = "") else "y"
  data <- if (isTRUE(data_overlay) && !isTRUE(deriv)) {
    data.frame(setNames(list(object$xz[, 1L]), payload$xname),
               y = object$y)
  } else {
    NULL
  }

  if (!identical(plot.behavior, "data"))
    .crs_plot_curve_draw(frame = frame, xname = payload$xname, yname = yname,
                         ci = ci, data = data,
                         data_overlay = data_overlay, ...)
  if (!identical(plot.behavior, "plot")) return(frame)
  invisible(frame)
}

.crs_plot_iv_deriv_public <- function(object,
                                      plot.call,
                                      plot.data = FALSE,
                                      phi = FALSE,
                                      ...) {
  raw.dots <- plot.call$...
  .crs_plot_validate_public_dots(raw.dots, context = "plot.crsivderiv",
                                 extra = c("plot.data", "phi"))
  dots <- list(...)
  .crs_plot_reject_curve_controls(dots, context = "plot.crsivderiv",
                                  allow.data.overlay = TRUE)
  plot.behavior <- .crs_plot_curve_behavior(dots)
  if (!is.null(dots$output)) dots$output <- NULL
  if (!is.null(dots$plot.behavior)) dots$plot.behavior <- NULL
  data_overlay <- if (!is.null(dots$data_overlay)) {
    .crs_plot_match_flag(dots$data_overlay, "data_overlay")
  } else {
    plot.data
  }
  dots$data_overlay <- NULL

  payload <- .crs_plot_payload_iv_deriv(object, phi = phi)
  frame <- payload$data
  names(frame)[1L] <- payload$xname
  yname <- if (isTRUE(phi)) "y" else paste("d", "y", "/d", payload$xname, sep = "")
  data <- if (isTRUE(data_overlay) && isTRUE(phi)) {
    data.frame(setNames(list(object$xz[, 1L]), payload$xname),
               y = object$y)
  } else {
    NULL
  }

  if (!identical(plot.behavior, "data"))
    .crs_plot_curve_draw(frame = frame, xname = payload$xname, yname = yname,
                         ci = FALSE, data = data,
                         data_overlay = isTRUE(data_overlay) && isTRUE(phi),
                         ...)
  if (!identical(plot.behavior, "plot")) return(frame)
  invisible(frame)
}

.crs_plot_clsd_public <- function(object,
                                  plot.call,
                                  er = TRUE,
                                  distribution = FALSE,
                                  derivative = FALSE,
                                  ...) {
  raw.dots <- plot.call$...
  .crs_plot_validate_public_dots(raw.dots, context = "plot.clsd",
                                 extra = c("er", "distribution", "derivative"))
  dots <- list(...)
  .crs_plot_reject_curve_controls(dots, context = "plot.clsd",
                                  allow.data.overlay = FALSE)
  plot.behavior <- .crs_plot_curve_behavior(dots)
  if (!is.null(dots$output)) dots$output <- NULL
  if (!is.null(dots$plot.behavior)) dots$plot.behavior <- NULL

  payload <- .crs_plot_payload_clsd(object, er = er,
                                    distribution = distribution,
                                    derivative = derivative)
  frame <- payload$data
  names(frame) <- c("x", "fit")
  yname <- switch(payload$view,
                  distribution = "Distribution",
                  derivative = "Density Derivative",
                  density = "Density",
                  "Fit")

  if (!identical(plot.behavior, "data"))
    .crs_plot_curve_draw(frame = frame, xname = "x", yname = yname,
                         ci = FALSE, data = NULL,
                         data_overlay = FALSE, ...)
  if (!identical(plot.behavior, "plot")) return(frame)
  invisible(frame)
}
