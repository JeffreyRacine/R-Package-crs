## Internal plot payload builders.
##
## These helpers are the first step of the CRS plot modernization path: they
## normalize object-fed plot data before any renderer owns public behavior.

.crs_plot_xq_vector <- function(object, xq) {
  p <- NCOL(object$xz)
  if (length(xq) == 1L) return(rep(xq, p))
  if (length(xq) != p) stop("xq must be scalar or have one value per predictor")
  xq
}

.crs_plot_slice_newdata <- function(object, slice.index, num.eval, xtrim, xq) {
  xz <- object$xz
  p <- NCOL(xz)
  if (slice.index < 1L || slice.index > p) stop("invalid plot slice index")

  if (!is.factor(xz[, slice.index])) {
    xlim <- trim.quantiles(xz[, slice.index], xtrim)
    slice.values <- seq(xlim[1L], xlim[2L], length = num.eval)
    neval <- num.eval
  } else {
    slice.values <- factor(levels(xz[, slice.index]),
                           levels = levels(xz[, slice.index]),
                           ordered = is.ordered(xz[, slice.index]))
    neval <- length(slice.values)
  }

  newdata <- vector("list", p)
  for (j in seq_len(p)) {
    if (j == slice.index) {
      newdata[[j]] <- slice.values
    } else if (!is.factor(xz[, j])) {
      newdata[[j]] <- rep(uocquantile(xz[, j], prob = xq[j]), neval)
    } else {
      newdata[[j]] <- factor(rep(uocquantile(xz[, j], prob = xq[j]), neval),
                             levels = levels(xz[, j]),
                             ordered = is.ordered(xz[, j]))
    }
  }

  newdata <- as.data.frame(newdata)
  names(newdata) <- names(xz)
  newdata
}

.crs_plot_surface_newdata <- function(object, num.eval, xtrim) {
  if (!is.null(object$num.z))
    stop("2D plot surfaces are supported only for continuous predictors")
  if (!identical(object$num.x, 2L))
    stop("2D plot surfaces are supported only for two continuous predictors")

  xz <- object$xz
  xlim <- trim.quantiles(xz[, 1L], xtrim)
  ylim <- trim.quantiles(xz[, 2L], xtrim)
  x1.seq <- seq(xlim[1L], xlim[2L], length = num.eval)
  x2.seq <- seq(ylim[1L], ylim[2L], length = num.eval)
  grid <- expand.grid(x1.seq, x2.seq)
  newdata <- data.frame(grid[, 1L], grid[, 2L])
  names(newdata) <- names(xz)
  list(newdata = newdata, x = x1.seq, y = x2.seq)
}

.crs_plot_prediction_frame <- function(object, newdata, deriv = 0L, ci = FALSE) {
  fit <- predict(object, newdata = newdata, deriv = deriv)
  out <- data.frame(newdata, fit = as.numeric(fit), check.names = FALSE)

  if (isTRUE(ci)) {
    lwr <- attr(fit, "lwr")
    upr <- attr(fit, "upr")
    if (!is.null(lwr) && !is.null(upr)) {
      out$lwr <- as.numeric(lwr)
      out$upr <- as.numeric(upr)
    }
  }

  out
}

.crs_plot_payload_regression <- function(object,
                                         deriv = 0L,
                                         gradients = FALSE,
                                         gradient.order = 1L,
                                         ci = FALSE,
                                         num.eval = 100L,
                                         neval = NULL,
                                         xtrim = 0,
                                         xq = 0.5,
                                         perspective = FALSE,
                                         legacy = FALSE,
                                         display.nomad.progress = FALSE,
                                         display.warnings = TRUE) {
  if (!inherits(object, "crs")) stop("object must inherit from class 'crs'")
  if (!is.null(neval)) num.eval <- neval
  if (isTRUE(gradients)) deriv <- gradient.order
  if (!is.numeric(deriv) || length(deriv) != 1L || is.na(deriv) || deriv < 0)
    stop("deriv must be a non-negative scalar")
  if (!is.numeric(num.eval) || length(num.eval) != 1L || is.na(num.eval) ||
      num.eval < 2L)
    stop("num.eval must be a scalar integer >= 2")
  if (!is.numeric(xtrim) || length(xtrim) != 1L || is.na(xtrim) ||
      xtrim < 0 || xtrim >= 0.5)
    stop("xtrim must be in [0, 0.5)")

  num.eval <- as.integer(num.eval)

  if (isTRUE(legacy)) {
    errors <- if (isTRUE(ci)) "asymptotic" else "none"
    if (deriv > 0) {
      payload <- plot(object, gradients = TRUE, gradient_order = deriv,
                      errors = errors, neval = num.eval,
                      xtrim = xtrim, xq = xq, output = "data",
                      display.nomad.progress = display.nomad.progress,
                      display.warnings = display.warnings)
    } else {
      payload <- plot(object, errors = errors, neval = num.eval,
                      xtrim = xtrim, xq = xq, output = "data",
                      display.nomad.progress = display.nomad.progress,
                      display.warnings = display.warnings)
    }
    if (is.list(payload) && is.null(names(payload)) &&
        length(payload) == NCOL(object$xz))
      names(payload) <- names(object$xz)
    return(structure(list(route = "crs",
                          view = if (deriv > 0) "derivative" else "fit",
                          source = "legacy",
                          deriv = deriv,
                          ci = isTRUE(ci),
                          tau = object$tau,
                          perspective = FALSE,
                          slices = payload),
                     class = "crs_plot_payload"))
  }

  if (deriv > 0) {
    if (isTRUE(perspective))
      stop("2D derivative plot payloads are not implemented yet")
    raw.slices <- .crs_plot_derivative_slices(
      object = object,
      deriv = deriv,
      ci = ci,
      num.eval = num.eval,
      xtrim = xtrim,
      xq = xq,
      plot.errors.type = "standard",
      display.warnings = display.warnings
    )
    slices <- Map(function(slice, nm) {
      names(slice)[2L] <- "fit"
      slice
    }, raw.slices, names(raw.slices))
    names(slices) <- names(raw.slices)
    return(structure(list(route = "crs",
                          view = "derivative",
                          source = "payload",
                          deriv = deriv,
                          ci = isTRUE(ci),
                          tau = object$tau,
                          perspective = FALSE,
                          slices = slices),
                     class = "crs_plot_payload"))
  }

  if (isTRUE(perspective)) {
    surface <- .crs_plot_surface_newdata(object, num.eval = num.eval,
                                         xtrim = xtrim)
    frame <- .crs_plot_prediction_frame(object, surface$newdata,
                                        deriv = 0L, ci = ci)
    z <- matrix(frame$fit, num.eval, num.eval)
    return(structure(list(route = "crs",
                          view = "surface",
                          source = "payload",
                          deriv = 0L,
                          ci = isTRUE(ci),
                          tau = object$tau,
                          perspective = TRUE,
                          x = surface$x,
                          y = surface$y,
                          z = z,
                          data = frame),
                     class = "crs_plot_payload"))
  }

  xq <- .crs_plot_xq_vector(object, xq)
  slices <- vector("list", NCOL(object$xz))
  names(slices) <- names(object$xz)
  for (i in seq_len(NCOL(object$xz))) {
    newdata <- .crs_plot_slice_newdata(object, i, num.eval, xtrim, xq)
    frame <- .crs_plot_prediction_frame(object, newdata,
                                        deriv = deriv, ci = ci)
    slices[[i]] <- frame
  }

  structure(list(route = "crs",
                 view = if (deriv > 0) "derivative" else "fit",
                 source = "payload",
                 deriv = deriv,
                 ci = isTRUE(ci),
                 tau = object$tau,
                 perspective = FALSE,
                 slices = slices),
            class = "crs_plot_payload")
}

.crs_plot_payload_iv <- function(object,
                                 deriv = FALSE,
                                 ci = FALSE,
                                 errors = "none",
                                 xtrim = 0) {
  if (!inherits(object, "crsiv")) stop("object must inherit from class 'crsiv'")
  if (object$num.x > 1 || !is.null(object$num.z))
    stop("only univariate z is supported")
  if (!is.logical(deriv) || length(deriv) != 1L || is.na(deriv))
    stop("deriv must be TRUE/FALSE")
  if (!is.logical(ci) || length(ci) != 1L || is.na(ci))
    stop("ci must be TRUE/FALSE")
  errors <- .crs_plot_scalar_match(errors,
                                   c("none", "bootstrap", "asymptotic"),
                                   "errors")
  ci <- isTRUE(ci) || !identical(errors, "none")
  if (!is.numeric(xtrim) || length(xtrim) != 1L || is.na(xtrim) ||
      xtrim < 0 || xtrim >= 0.5)
    stop("xtrim must be in [0, 0.5)")

  z <- object$xz[, 1L]
  xlim <- if (xtrim == 0) {
    range(z, na.rm = TRUE)
  } else {
    as.numeric(stats::quantile(z, probs = c(xtrim, 1 - xtrim),
                               names = FALSE, na.rm = TRUE, type = 7))
  }
  keep <- (z >= xlim[1L]) & (z <= xlim[2L])
  ord <- order(z[keep])
  z.plot <- z[keep][ord]

  if (isTRUE(deriv)) {
    if (is.null(object$deriv.mat))
      stop("deriv.mat not found: was crsiv called with deriv > 0?")
    fit <- object$deriv.mat[keep, 1L][ord]
    lwr <- if (isTRUE(ci) && !is.null(object$deriv.mat.lwr))
      object$deriv.mat.lwr[keep, 1L][ord] else NULL
    upr <- if (isTRUE(ci) && !is.null(object$deriv.mat.upr))
      object$deriv.mat.upr[keep, 1L][ord] else NULL
  } else {
    fit <- object$phi[keep][ord]
    lwr <- if (isTRUE(ci) && !is.null(object$phi.lwr))
      object$phi.lwr[keep][ord] else NULL
    upr <- if (isTRUE(ci) && !is.null(object$phi.upr))
      object$phi.upr[keep][ord] else NULL
  }

  frame <- data.frame(z = z.plot, fit = fit)
  if (!is.null(lwr) && !is.null(upr)) {
    frame$lwr <- lwr
    frame$upr <- upr
  }

  structure(list(route = "crsiv",
                 view = if (isTRUE(deriv)) "derivative" else "fit",
                 ci = isTRUE(ci),
                 data = frame,
                 xlim = xlim,
                 xname = object$xnames[1L]),
            class = "crs_plot_payload")
}

.crs_plot_payload_iv_deriv <- function(object, phi = FALSE) {
  if (!inherits(object, "crsivderiv"))
    stop("object must inherit from class 'crsivderiv'")
  if (object$num.x > 1 || !is.null(object$num.z))
    stop("only univariate z is supported")
  if (!is.logical(phi) || length(phi) != 1L || is.na(phi))
    stop("phi must be TRUE/FALSE")

  z <- object$xz[, 1L]
  fit <- if (isTRUE(phi)) object$phi else object$phi.prime
  ord <- order(z)

  structure(list(route = "crsivderiv",
                 view = if (isTRUE(phi)) "fit" else "derivative",
                 data = data.frame(z = z[ord], fit = fit[ord]),
                 xname = object$xnames[1L]),
            class = "crs_plot_payload")
}

.crs_plot_payload_clsd <- function(object,
                                   er = TRUE,
                                   distribution = FALSE,
                                   derivative = FALSE) {
  if (!inherits(object, "clsd")) stop("object must inherit from class 'clsd'")
  if (!is.logical(er) || length(er) != 1L || is.na(er))
    stop("er must be TRUE/FALSE")
  if (!is.logical(distribution) || length(distribution) != 1L ||
      is.na(distribution))
    stop("distribution must be TRUE/FALSE")
  if (!is.logical(derivative) || length(derivative) != 1L ||
      is.na(derivative))
    stop("derivative must be TRUE/FALSE")

  if (isTRUE(er)) {
    order.x <- order(object$xer)
    x <- object$xer[order.x]
    if (isTRUE(distribution)) {
      y <- object$distribution.er[order.x]
      view <- "distribution"
    } else if (isTRUE(derivative)) {
      y <- object$density.deriv.er[order.x]
      view <- "derivative"
    } else {
      y <- object$density.er[order.x]
      view <- "density"
    }
  } else {
    order.x <- order(object$x)
    x <- object$x[order.x]
    if (isTRUE(distribution)) {
      y <- object$distribution[order.x]
      view <- "distribution"
    } else if (isTRUE(derivative)) {
      y <- object$density.deriv[order.x]
      view <- "derivative"
    } else {
      y <- object$density[order.x]
      view <- "density"
    }
  }

  structure(list(route = "clsd",
                 view = view,
                 er = isTRUE(er),
                 data = data.frame(x = x, y = y)),
            class = "crs_plot_payload")
}
