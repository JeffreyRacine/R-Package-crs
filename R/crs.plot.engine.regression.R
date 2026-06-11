## Internal CRS regression plot engines.

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

.crs_plot_payload_to_derivative_slices <- function(payload) {
  out <- .crs_plot_payload_to_legacy_slices(payload)
  for (i in seq_along(out)) {
    names(out[[i]])[2L] <- "deriv"
  }
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
                                           deriv = 0L,
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
    ylab <- if (deriv > 0L) {
      if (!is.factor(object$xz[[nm]])) {
        paste("Order", deriv, "Derivative")
      } else {
        "Difference in Levels"
      }
    } else if (is.null(object$tau)) {
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

    if (isTRUE(data_overlay) && identical(as.integer(deriv), 0L) &&
        nm %in% names(object$xz)) {
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

.crs_plot_mean_bootstrap_slices <- function(object,
                                            num.eval,
                                            xtrim,
                                            xq,
                                            plot.errors.boot.num,
                                            plot.errors.type,
                                            plot.errors.alpha,
                                            display.nomad.progress,
                                            display.warnings) {
  xq <- .crs_plot_xq_vector(object, xq)
  slices <- vector("list", NCOL(object$xz))
  names(slices) <- names(object$xz)

  for (i in seq_len(NCOL(object$xz))) {
    newdata <- .crs_plot_slice_newdata(object, i, num.eval, xtrim, xq)
    target.label <- .crs_plot_regression_bootstrap_target_label(
      object = object,
      slice.index = i,
      gradients = FALSE
    )
    boot <- .crs.bootstrap.matrix(object = object,
                                  newdata = newdata,
                                  deriv = 0L,
                                  deriv.index = i,
                                  boot.num = plot.errors.boot.num,
                                  display.warnings = display.warnings,
                                  display.nomad.progress = display.nomad.progress,
                                  progress.target = target.label)
    if (identical(plot.errors.type, "all")) {
      all.bounds <- .crs.bootstrap.bounds(boot$boot.mat,
                                          plot.errors.alpha,
                                          "all",
                                          boot$center)
      slices[[i]] <- data.frame(newdata[, i],
                                boot$center,
                                all.bounds$pointwise[, 1L],
                                all.bounds$pointwise[, 2L],
                                all.bounds$simultaneous[, 1L],
                                all.bounds$simultaneous[, 2L],
                                all.bounds$bonferroni[, 1L],
                                all.bounds$bonferroni[, 2L])
      names(slices[[i]]) <- c(names(newdata)[i], "mean", "lwr", "upr",
                              "lwr.sim", "upr.sim", "lwr.bonf", "upr.bonf")
    } else {
      bounds <- .crs.bootstrap.bounds(boot$boot.mat,
                                      plot.errors.alpha,
                                      plot.errors.type,
                                      boot$center)
      slices[[i]] <- data.frame(newdata[, i], boot$center, bounds)
      names(slices[[i]]) <- c(names(newdata)[i], "mean", "lwr", "upr")
    }
  }

  slices
}

.crs_plot_derivative_slices <- function(object,
                                        deriv,
                                        ci,
                                        num.eval,
                                        xtrim,
                                        xq,
                                        plot.errors.type,
                                        display.warnings = TRUE) {
  if (!inherits(object, "crs")) stop("object must inherit from class 'crs'")
  if (!is.numeric(deriv) || length(deriv) != 1L || is.na(deriv) || deriv <= 0)
    stop("deriv must be a positive scalar for derivative plot slices")

  basis <- object$basis
  prune <- object$prune
  prune.index <- object$prune.index
  xz <- object$xz
  y <- object$y

  if (!object$kernel) {
    xztmp <- splitFrame(xz)
  } else {
    xztmp <- splitFrame(xz, factor.to.numeric = TRUE)
  }
  x <- xztmp$x
  z <- xztmp$z
  is.ordered.z <- xztmp$is.ordered.z

  knots <- object$knots
  K <- object$K
  degree <- object$degree
  include <- object$include
  lambda <- object$lambda
  tau <- object$tau
  weights <- object$weights
  xq <- .crs_plot_xq_vector(object, xq)

  slices <- vector("list", NCOL(object$xz))
  names(slices) <- names(object$xz)
  m <- 0L
  i.numeric <- 0L

  for (i in seq_len(NCOL(object$xz))) {
    if (!is.factor(object$xz[, i])) {
      i.numeric <- i.numeric + 1L
      newdata <- matrix(NA, nrow = num.eval, ncol = NCOL(object$xz))
      neval <- num.eval
      m <- m + 1L
    } else {
      newdata <- matrix(NA,
                        nrow = length(levels(object$xz[, i])),
                        ncol = NCOL(object$xz))
      neval <- length(levels(object$xz[, i]))
    }

    newdata <- data.frame(newdata)
    newdata.base <- data.frame(newdata)

    if (!is.factor(object$xz[, i])) {
      xlim <- trim.quantiles(object$xz[, i], xtrim)
      newdata[, i] <- seq(xlim[1L], xlim[2L], length = neval)
    } else {
      newdata[, i] <- factor(levels(object$xz[, i]),
                             levels = levels(object$xz[, i]),
                             ordered = is.ordered(object$xz[, i]))
      newdata.base[, i] <- factor(rep(levels(object$xz[, i])[1L], neval),
                                  levels = levels(object$xz[, i]),
                                  ordered = is.ordered(object$xz[, i]))
    }

    for (j in (seq_len(NCOL(object$xz)))[-i]) {
      if (!is.factor(object$xz[, j])) {
        newdata[, j] <- rep(uocquantile(object$xz[, j], prob = xq[j]), neval)
        newdata.base[, j] <- rep(uocquantile(object$xz[, j], prob = xq[j]),
                                 neval)
      } else {
        newdata[, j] <- factor(rep(uocquantile(object$xz[, j], prob = xq[j]),
                                   neval),
                               levels = levels(object$xz[, j]),
                               ordered = is.ordered(object$xz[, j]))
        newdata.base[, j] <- factor(rep(uocquantile(object$xz[, j],
                                                    prob = xq[j]), neval),
                                    levels = levels(object$xz[, j]),
                                    ordered = is.ordered(object$xz[, j]))
      }
    }

    newdata <- data.frame(newdata)
    names(newdata) <- names(object$xz)
    newdata.base <- data.frame(newdata.base)
    names(newdata.base) <- names(object$xz)

    if (!object$kernel) {
      xztmp <- splitFrame(data.frame(newdata))
    } else {
      xztmp <- splitFrame(data.frame(newdata), factor.to.numeric = TRUE)
    }
    xeval <- xztmp$x
    zeval <- xztmp$z
    is.ordered.z <- xztmp$is.ordered.z

    if (!object$kernel) {
      xztmp <- splitFrame(data.frame(newdata.base))
    } else {
      xztmp <- splitFrame(data.frame(newdata.base), factor.to.numeric = TRUE)
    }
    xeval.base <- xztmp$x
    zeval.base <- xztmp$z
    is.ordered.z <- xztmp$is.ordered.z

    if (!object$kernel) {
      if (!is.factor(newdata[, i])) {
        if (deriv <= degree[i.numeric]) {
          tmp <- derivFactorSpline(x = x,
                                   y = y,
                                   z = z,
                                   K = K,
                                   I = include,
                                   xeval = xeval,
                                   zeval = zeval,
                                   knots = knots,
                                   basis = basis,
                                   deriv.index = m,
                                   deriv = deriv,
                                   prune.index = prune.index,
                                   tau = tau,
                                   weights = weights)
        } else {
          tmp <- matrix(0, nrow(newdata), 3L)
        }
        deriv.est <- tmp[, 1L]
        deriv.lwr <- tmp[, 2L]
        deriv.upr <- tmp[, 3L]
      } else {
        zpred <- preditFactorSpline(x = x, y = y, z = z, K = K, I = include,
                                    xeval = xeval, zeval = zeval,
                                    knots = knots, basis = basis,
                                    prune = prune,
                                    prune.index = prune.index,
                                    tau = tau, weights = weights)$fitted.values
        zpred.base <- preditFactorSpline(x = x, y = y, z = z, K = K,
                                         I = include,
                                         xeval = xeval.base,
                                         zeval = zeval.base,
                                         knots = knots, basis = basis,
                                         prune = prune,
                                         prune.index = prune.index,
                                         tau = tau,
                                         weights = weights)$fitted.values
        deriv.est <- zpred[, 1L] - zpred.base[, 1L]
        deriv.lwr <- deriv.est -
          stats::qnorm(0.975) * sqrt(zpred[, 4L]^2 + zpred.base[, 4L]^2)
        deriv.upr <- deriv.est +
          stats::qnorm(0.975) * sqrt(zpred[, 4L]^2 + zpred.base[, 4L]^2)
      }
    } else {
      if (!is.factor(newdata[, i])) {
        if (deriv <= degree[i.numeric]) {
          tmp <- derivKernelSpline(x = x,
                                   y = y,
                                   z = z,
                                   K = K,
                                   lambda = lambda,
                                   is.ordered.z = is.ordered.z,
                                   xeval = xeval,
                                   zeval = zeval,
                                   knots = knots,
                                   basis = basis,
                                   deriv.index = m,
                                   deriv = deriv,
                                   tau = tau,
                                   weights = weights)
        } else {
          tmp <- matrix(0, nrow(newdata), 3L)
        }
        deriv.est <- tmp[, 1L]
        deriv.lwr <- tmp[, 2L]
        deriv.upr <- tmp[, 3L]
      } else {
        z <- as.matrix(z)
        zeval <- as.matrix(zeval)
        zeval.base <- as.matrix(zeval.base)
        zpred <- predictKernelSpline(x = x, y = y, z = z, K = K,
                                     lambda = lambda,
                                     is.ordered.z = is.ordered.z,
                                     xeval = xeval, zeval = zeval,
                                     knots = knots, basis = basis,
                                     tau = tau,
                                     weights = weights)$fitted.values
        zpred.base <- predictKernelSpline(x = x, y = y, z = z, K = K,
                                          lambda = lambda,
                                          is.ordered.z = is.ordered.z,
                                          xeval = xeval.base,
                                          zeval = zeval.base,
                                          knots = knots, basis = basis,
                                          tau = tau,
                                          weights = weights)$fitted.values
        deriv.est <- zpred[, 1L] - zpred.base[, 1L]
        deriv.lwr <- deriv.est -
          stats::qnorm(0.975) * sqrt(zpred[, 4L]^2 + zpred.base[, 4L]^2)
        deriv.upr <- deriv.est +
          stats::qnorm(0.975) * sqrt(zpred[, 4L]^2 + zpred.base[, 4L]^2)
      }
    }

    if (!isTRUE(ci)) {
      slices[[i]] <- data.frame(newdata[, i], deriv.est)
      names(slices[[i]]) <- c(names(newdata)[i], "deriv")
    } else if (identical(plot.errors.type, "all")) {
      slices[[i]] <- data.frame(newdata[, i], deriv.est,
                                deriv.lwr, deriv.upr,
                                deriv.lwr, deriv.upr,
                                deriv.lwr, deriv.upr)
      names(slices[[i]]) <- c(names(newdata)[i], "deriv", "lwr", "upr",
                              "lwr.sim", "upr.sim", "lwr.bonf", "upr.bonf")
      if (is.factor(newdata[, i]) && isTRUE(display.warnings))
        warning("bootstrap-all for factor derivatives currently reuses standard bounds for this slice")
    } else {
      slices[[i]] <- data.frame(newdata[, i], deriv.est, deriv.lwr, deriv.upr)
      names(slices[[i]]) <- c(names(newdata)[i], "deriv", "lwr", "upr")
    }
  }

  slices
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
  deriv <- as.numeric(.crs_plot_scalar_default(dots$deriv, 0L))
  num.eval <- as.integer(.crs_plot_scalar_default(dots$num.eval, 100L))
  xtrim <- .crs_plot_scalar_default(dots$xtrim, 0)
  xq <- .crs_plot_scalar_default(dots$xq, 0.5)
  common.scale <- isTRUE(.crs_plot_scalar_default(dots$common.scale, TRUE))
  data_overlay <- isTRUE(.crs_plot_scalar_default(dots$plot.data.overlay, TRUE))
  perspective <- isTRUE(.crs_plot_scalar_default(dots$perspective, FALSE))
  if (isTRUE(perspective))
    stop("modern 2D regression plot route is not implemented yet",
         call. = FALSE)
  plot.errors.method <- .crs_plot_scalar_default(dots$plot.errors.method,
                                                 "none")
  ci <- isTRUE(.crs_plot_scalar_default(
    dots$ci, !identical(plot.errors.method, "none")
  ))
  plot.errors.type <- .crs_plot_scalar_default(dots$plot.errors.type,
                                               "standard")
  plot.errors.alpha <- .crs_plot_scalar_default(dots$plot.errors.alpha, 0.05)
  plot.errors.boot.num <- as.integer(.crs_plot_scalar_default(
    dots$plot.errors.boot.num, 99L
  ))

  if (deriv > 0) {
    if (isTRUE(ci) && identical(plot.errors.method, "bootstrap"))
      stop("bootstrap intervals for derivative plots are not implemented",
           call. = FALSE)
    slices <- .crs_plot_derivative_slices(
      object = object,
      deriv = deriv,
      ci = ci,
      num.eval = num.eval,
      xtrim = xtrim,
      xq = xq,
      plot.errors.type = plot.errors.type,
      display.warnings = .crs_plot_scalar_default(dots$display.warnings, TRUE)
    )
  } else if (isTRUE(ci) && identical(plot.errors.method, "bootstrap")) {
    slices <- .crs_plot_mean_bootstrap_slices(
      object = object,
      num.eval = num.eval,
      xtrim = xtrim,
      xq = xq,
      plot.errors.boot.num = plot.errors.boot.num,
      plot.errors.type = plot.errors.type,
      plot.errors.alpha = plot.errors.alpha,
      display.nomad.progress = .crs_plot_scalar_default(
        dots$display.nomad.progress, FALSE
      ),
      display.warnings = .crs_plot_scalar_default(dots$display.warnings, TRUE)
    )
  } else {
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
  }

  if (!identical(plot.behavior, "data")) {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar), add = TRUE)
    render.dots <- dots[setdiff(names(dots),
                                c("plot.behavior", "plot.data.overlay",
                                  "num.eval", "xtrim", "xq", "ci",
                                  "common.scale", "deriv",
                                  "plot.errors.method", "plot.errors.type",
                                  "plot.errors.alpha",
                                  "plot.errors.boot.num",
                                  "display.nomad.progress",
                                  "display.warnings"))]
    do.call(.crs_plot_render_regression_1d,
            c(list(object = object,
                   slices = slices,
                   deriv = deriv,
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
                                           ...) {
  if (!inherits(object, "crs")) stop("object must inherit from class 'crs'")

  raw.dots <- plot.call$...
  dot.names <- names(raw.dots)
  if (is.null(dot.names)) dot.names <- character()
  .crs_plot_validate_public_dots(raw.dots, context = "plot.crs")
  dots <- list(...)
  dots <- .crs_plot_normalize_public_dots(dots, context = "plot.crs")

  plot.behavior <- if (!is.null(dots$plot.behavior)) {
    match.arg(dots$plot.behavior, c("plot", "plot-data", "data"))
  } else {
    "plot"
  }
  gradients <- isTRUE(.crs_plot_scalar_default(dots$gradients, FALSE))
  gradient.order <- .crs_plot_scalar_default(dots$gradient.order, 1L)
  if (!is.numeric(gradient.order) || any(is.na(gradient.order)) ||
      any(gradient.order < 1L))
    stop("gradient_order must contain positive numeric values",
         call. = FALSE)
  if (length(gradient.order) != 1L)
    stop("plot.crs gradients currently require scalar gradient_order",
         call. = FALSE)
  deriv <- if (isTRUE(gradients)) as.integer(gradient.order) else 0L

  num.eval <- as.integer(.crs_plot_scalar_default(dots$num.eval, 50L))
  xtrim <- .crs_plot_scalar_default(dots$xtrim, 0)
  xq <- .crs_plot_scalar_default(dots$xq, 0.5)
  common.scale <- isTRUE(.crs_plot_scalar_default(dots$common.scale, TRUE))
  display.nomad.progress <- isTRUE(.crs_plot_scalar_default(
    dots$display.nomad.progress, TRUE
  ))
  display.warnings <- isTRUE(.crs_plot_scalar_default(
    dots$display.warnings, TRUE
  ))
  plot.errors.method <- .crs_plot_scalar_default(dots$plot.errors.method,
                                                 "none")
  plot.errors.method <- .crs_plot_scalar_match(plot.errors.method,
                                               c("none", "bootstrap",
                                                 "asymptotic"),
                                               "errors")
  plot.errors.type <- .crs_plot_scalar_default(dots$plot.errors.type,
                                               "standard")
  plot.errors.alpha <- .crs_plot_scalar_default(dots$plot.errors.alpha, 0.05)
  plot.errors.boot.num <- as.integer(.crs_plot_scalar_default(
    dots$plot.errors.boot.num, 99L
  ))
  plot.errors.boot.method <- .crs_plot_scalar_default(
    dots$plot.errors.boot.method, "inid"
  )
  plot.errors.center <- .crs_plot_scalar_default(dots$plot.errors.center,
                                                 "estimate")
  if (!identical(plot.errors.center, "estimate"))
    stop("plot.crs currently supports center=\"estimate\" only",
         call. = FALSE)
  if (identical(plot.errors.method, "bootstrap") &&
      !identical(plot.errors.boot.method, "inid"))
    stop("plot.crs bootstrap intervals currently support bootstrap=\"inid\" only",
         call. = FALSE)
  if (identical(plot.errors.method, "asymptotic") &&
      !identical(plot.errors.type, "standard"))
    stop("plot.crs asymptotic intervals currently support band=\"pmzsd\" only",
         call. = FALSE)

  ci <- !identical(plot.errors.method, "none")
  perspective <- isTRUE(.crs_plot_scalar_default(dots$perspective, TRUE))
  renderer <- .crs_plot_scalar_default(dots$renderer, "base")
  renderer <- match.arg(renderer, c("base", "rgl"))
  surface.supported <- is.null(object$num.z) && identical(object$num.x, 2L)
  if ("perspective" %in% dot.names && isTRUE(perspective) &&
      !isTRUE(surface.supported) && !isTRUE(gradients))
    stop("2D plot surfaces are supported only for two continuous predictors",
         call. = FALSE)
  surface.request <- isTRUE(perspective) && isTRUE(surface.supported) &&
    !isTRUE(gradients)
  if ("renderer" %in% dot.names && !isTRUE(surface.request))
    stop("renderer is supported only for 2D fitted-function surfaces",
         call. = FALSE)

  bridge <- dots
  bridge$plot.behavior <- plot.behavior
  bridge$ci <- ci
  bridge$deriv <- deriv
  bridge$num.eval <- num.eval
  bridge$xtrim <- xtrim
  bridge$xq <- xq
  bridge$common.scale <- common.scale
  bridge$display.nomad.progress <- display.nomad.progress
  bridge$display.warnings <- display.warnings
  bridge$plot.errors.method <- plot.errors.method
  bridge$plot.errors.type <- plot.errors.type
  bridge$plot.errors.alpha <- plot.errors.alpha
  if (identical(plot.errors.method, "bootstrap"))
    bridge$plot.errors.boot.num <- plot.errors.boot.num

  if (isTRUE(ci) && isTRUE(gradients) &&
      identical(plot.errors.method, "bootstrap"))
    stop("bootstrap intervals for derivative plots are not implemented",
         call. = FALSE)

  if (isTRUE(surface.request)) {
    if (isTRUE(ci))
      stop("plot.crs surface route currently supports point estimates only",
           call. = FALSE)
    bridge$renderer <- renderer
    return(do.call(.crs_plot_regression_surface_shadow,
                   c(list(object = object, .plot_dots_call = raw.dots),
                     bridge)))
  }

  do.call(.crs_plot_regression_1d_shadow,
          c(list(object = object, .plot_dots_call = raw.dots), bridge))
}
