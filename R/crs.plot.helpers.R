## CRS plot helpers copied/adapted from the modern np plot layer.
##
## The public plot argument grammar is intentionally NP-shaped. CRS-specific
## internals translate those public controls into spline prediction payloads.

.crs_plot_scalar_match <- function(value, choices, argname) {
  if (is.null(value)) return(NULL)
  if (length(value) != 1L || is.na(value))
    stop(sprintf("%s must be one of %s",
                 argname,
                 paste(sprintf("\"%s\"", choices), collapse = ", ")),
         call. = FALSE)
  value <- as.character(value)
  if (!(value %in% choices))
    stop(sprintf("%s must be one of %s",
                 argname,
                 paste(sprintf("\"%s\"", choices), collapse = ", ")),
         call. = FALSE)
  value
}

np_boot_control <- function(nonfixed = c("exact", "frozen"),
                            wild = c("rademacher", "mammen"),
                            blocklen = NULL) {
  nonfixed <- match.arg(nonfixed)
  wild <- match.arg(wild)
  if (!is.null(blocklen) &&
      (!is.numeric(blocklen) || length(blocklen) != 1L ||
       is.na(blocklen) || blocklen <= 0))
    stop("blocklen must be a positive numeric scalar", call. = FALSE)
  structure(list(nonfixed = nonfixed, wild = wild, blocklen = blocklen),
            class = "np_boot_control")
}

np_grid_control <- function(xtrim = NULL, xq = NULL, slices = NULL) {
  if (!is.null(xtrim) &&
      (!is.numeric(xtrim) || length(xtrim) != 2L ||
       any(is.na(xtrim)) || any(xtrim < 0) || any(xtrim > 1) ||
       xtrim[1L] >= xtrim[2L]))
    stop("xtrim must be a numeric length-two vector with 0 <= xtrim[1] < xtrim[2] <= 1",
         call. = FALSE)
  structure(list(xtrim = xtrim, xq = xq, slices = slices),
            class = "np_grid_control")
}

np_render_control <- function(style = c("band", "bar"),
                              bar = c("|", "I"),
                              bar_num = NULL) {
  style <- match.arg(style)
  bar <- match.arg(bar)
  if (!is.null(bar_num) &&
      (!is.numeric(bar_num) || length(bar_num) != 1L ||
       is.na(bar_num) || bar_num < 1))
    stop("bar_num must be a positive numeric scalar", call. = FALSE)
  structure(list(style = style, bar = bar, bar_num = bar_num),
            class = "np_render_control")
}

crs_boot_control <- function(...) {
  out <- np_boot_control(...)
  class(out) <- c("crs_boot_control", class(out))
  out
}

crs_grid_control <- function(...) {
  out <- np_grid_control(...)
  class(out) <- c("crs_grid_control", class(out))
  out
}

crs_render_control <- function(...) {
  out <- np_render_control(...)
  class(out) <- c("crs_render_control", class(out))
  out
}

.crs_boot_control <- crs_boot_control
.crs_grid_control <- crs_grid_control
.crs_render_control <- crs_render_control

.crs_plot_dot_names <- function(dots_call) {
  if (is.null(dots_call) || length(dots_call) == 0L)
    return(character())
  nms <- names(dots_call)
  if (is.null(nms))
    rep.int("", length(dots_call))
  else
    nms
}

.crs_plot_match_flag <- function(value, argname) {
  if (is.null(value)) return(NULL)
  if (!is.logical(value) || length(value) != 1L || is.na(value))
    stop(sprintf("%s must be TRUE or FALSE", argname), call. = FALSE)
  value
}

.crs_plot_scalar_default <- function(value, default) {
  if (is.null(value)) default else value
}

.crs_plot_viridis_at <- function(x, begin = 0.12, end = 0.88, alpha = 1) {
  x.names <- names(x)
  x <- as.double(x)
  x <- pmax.int(0, pmin.int(1, x))
  grid <- grDevices::hcl.colors(256L, palette = "viridis")
  pos <- begin + x * (end - begin)
  idx <- 1L + round(pos * (length(grid) - 1L))
  idx <- pmax.int(1L, pmin.int(length(grid), idx))
  cols <- grid[idx]
  if (!identical(alpha, 1))
    cols <- grDevices::adjustcolor(cols, alpha.f = alpha)
  names(cols) <- x.names
  cols
}

.crs_plot_viridis_role <- function(x) {
  unname(.crs_plot_viridis_at(x))[1L]
}

.crs_plot_color <- function(role, alpha = NULL) {
  spec <- switch(
    as.character(role)[1L],
    primary = list(col = .crs_plot_viridis_role(0.02), alpha = 1),
    median = list(col = .crs_plot_viridis_role(0.02), alpha = 1),
    secondary = list(col = .crs_plot_viridis_role(0.78), alpha = 1),
    fit = list(col = .crs_plot_viridis_role(0.02), alpha = 1),
    interval = list(col = .crs_plot_viridis_role(0.78), alpha = 1),
    interval_context = list(col = .crs_plot_viridis_role(0.78), alpha = 1),
    data_overlay = list(col = .crs_plot_viridis_role(0.08), alpha = 0.35),
    support = list(col = .crs_plot_viridis_role(0.12), alpha = 0.60),
    support_floor = list(col = .crs_plot_viridis_role(0.18), alpha = 0.55),
    support_grid = list(col = .crs_plot_viridis_role(0.50), alpha = 0.45),
    component_context = list(col = .crs_plot_viridis_role(0.68), alpha = 1),
    surface_border = list(col = .crs_plot_viridis_role(0.02), alpha = 1),
    context_wire = list(col = .crs_plot_viridis_role(0.25), alpha = 1),
    context_border = list(col = .crs_plot_viridis_role(0.30), alpha = 1),
    legend_bg = list(col = .crs_plot_viridis_role(0.98), alpha = 0.18),
    list(col = as.character(role)[1L], alpha = 1)
  )
  alpha <- if (is.null(alpha)) spec$alpha else alpha
  if (isTRUE(all.equal(alpha, 1))) {
    spec$col
  } else {
    grDevices::adjustcolor(spec$col, alpha.f = alpha)
  }
}

.crs_plot_pch <- function(role) {
  switch(as.character(role)[1L],
         data_overlay = 20L,
         fit = 1L,
         20L)
}

.crs_plot_cex <- function(role) {
  switch(as.character(role)[1L],
         data_overlay = 0.5,
         legend = 0.8,
         1)
}

.crs_plot_lwd <- function(role, base = graphics::par()$lwd) {
  base <- as.double(base)[1L]
  switch(as.character(role)[1L],
         primary = base,
         surface_border = 0.8 * base,
         interval = base,
         band_all_1d = 2 * base,
         band_all_surface = 2.15 * base,
         interval_surface = 2 * base,
         support = 1.25 * base,
         support_floor = 2 * base,
         support_grid = 0.9 * base,
         quantile_multi = 1.5 * base,
         component_context = base,
         base)
}

.crs_plot_lty <- function(role) {
  switch(as.character(role)[1L],
         solid = 1L,
         interval = 2L,
         center = 3L,
         lower_quantile = 2L,
         upper_quantile = 4L,
         pointwise = 2L,
         simultaneous = 3L,
         bonferroni = 4L,
         1L)
}

.crs_plot_surface_colors <- function(z, num.colors = 100L, alpha = 0.90) {
  z <- as.numeric(z)
  pal <- grDevices::hcl.colors(num.colors, palette = "viridis")
  if (!length(z) || all(!is.finite(z)))
    return(grDevices::adjustcolor(pal[1L], alpha.f = alpha))
  rng <- range(z, finite = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0)
    idx <- rep.int(ceiling(num.colors / 2), length(z))
  else
    idx <- pmax(1L, pmin(num.colors,
                         as.integer(cut(z, breaks = num.colors,
                                        labels = FALSE))))
  grDevices::adjustcolor(pal[idx], alpha.f = alpha)
}

.crs_plot_persp_surface_colors <- function(z, col = NULL,
                                           num.colors = 1000L) {
  if (!is.null(col)) return(col)
  z <- as.matrix(z)
  z.range <- range(z, finite = TRUE)
  palette_fun <- function(n) grDevices::hcl.colors(as.integer(n),
                                                   palette = "viridis")
  if (!all(is.finite(z.range)))
    return(palette_fun(1L))
  if (nrow(z) < 2L || ncol(z) < 2L)
    return(palette_fun(1L))
  if (isTRUE(all.equal(z.range[1L], z.range[2L])))
    return(palette_fun(1L))
  z.facet <- (z[-1L, -1L, drop = FALSE] +
                z[-nrow(z), -1L, drop = FALSE] +
                z[-1L, -ncol(z), drop = FALSE] +
                z[-nrow(z), -ncol(z), drop = FALSE]) / 4
  colorlut <- palette_fun(num.colors)
  scaled <- 1L + floor((length(colorlut) - 1L) *
                         (z.facet - z.range[1L]) / diff(z.range))
  scaled[!is.finite(scaled)] <- 1L
  scaled <- pmax.int(1L, pmin.int(length(colorlut), scaled))
  as.vector(matrix(colorlut[scaled],
                   nrow = nrow(z.facet),
                   ncol = ncol(z.facet)))
}

.crs_plot_rgl_surface_colors <- function(z, col = NULL, num.colors = 1000L) {
  if (!is.null(col)) return(col)
  z <- as.matrix(z)
  colorlut <- grDevices::hcl.colors(as.integer(num.colors),
                                    palette = "viridis")
  finite.z <- z[is.finite(z)]
  if (!length(finite.z))
    return(matrix(colorlut[1L], nrow = nrow(z), ncol = ncol(z)))
  zmin <- min(finite.z)
  zmax <- max(finite.z)
  if (!is.finite(zmin) || !is.finite(zmax) || zmax <= zmin)
    return(matrix(colorlut[ceiling(num.colors / 2)],
                  nrow = nrow(z), ncol = ncol(z)))
  idx <- floor((num.colors - 1L) * (z - zmin) / (zmax - zmin)) + 1L
  idx[!is.finite(idx)] <- 1L
  idx <- pmax.int(1L, pmin.int(num.colors, idx))
  matrix(colorlut[idx], nrow = nrow(z), ncol = ncol(z))
}

.crs_plot_user_args <- function(dots, type) {
  if (is.null(dots) || !length(dots)) return(list())
  prefix <- paste0(type, ".")
  nms <- names(dots)
  if (is.null(nms)) return(list())

  direct.names <- switch(type,
    plot = unique(c(
      setdiff(names(formals(graphics::plot.default)), c("x", "y", "...")),
      names(graphics::par(no.readonly = TRUE)),
      "type", "lty", "lwd", "col", "pch", "cex", "main", "sub",
      "xlab", "ylab", "xlim", "ylim"
    )),
    persp = {
      persp.default <- utils::getS3method("persp", "default", optional = TRUE)
      if (is.null(persp.default)) {
        character()
      } else {
        setdiff(names(formals(persp.default)), c("x", "y", "z", "..."))
      }
    },
    character()
  )

  direct <- dots[intersect(nms, direct.names)]
  grouped <- dots[[type]]
  if (!is.null(grouped)) {
    if (!is.list(grouped))
      stop(sprintf("%s must be a list when supplied as a grouped plot argument",
                   type),
           call. = FALSE)
    direct <- c(direct, grouped)
  }

  prefixed <- dots[startsWith(nms, prefix)]
  if (length(prefixed)) {
    names(prefixed) <- substring(names(prefixed), nchar(prefix) + 1L)
    direct <- c(direct, prefixed)
  }

  direct
}

.crs_plot_merge_user_args <- function(base.args, user.args) {
  if (!length(user.args)) return(base.args)
  for (nm in names(user.args)) base.args[[nm]] <- user.args[[nm]]
  base.args
}

.crs_plot_merge_override_args <- function(base.args, override.args) {
  if (is.null(override.args) || !length(override.args))
    return(base.args)
  if (is.null(base.args) || !length(base.args))
    return(override.args)

  base.names <- names(base.args)
  override.names <- names(override.args)
  dup <- intersect(override.names[!(is.na(override.names) |
                                      override.names == "")],
                   base.names[!(is.na(base.names) |
                                  base.names == "")])
  if (length(dup)) {
    keep <- is.na(base.names) | base.names == "" | !(base.names %in% dup)
    base.args <- base.args[keep]
  }
  c(base.args, override.args)
}

.crs_plot_overlay_range <- function(existing.range, y) {
  if (is.null(y)) return(existing.range)
  yr <- range(y, finite = TRUE)
  if (!length(yr) || any(!is.finite(yr))) return(existing.range)
  if (is.null(existing.range) || length(existing.range) < 2L ||
      all(!is.finite(existing.range))) {
    return(yr)
  }
  c(min(existing.range[1L], yr[1L], na.rm = TRUE),
    max(existing.range[2L], yr[2L], na.rm = TRUE))
}

.crs_plot_overlay_points_1d <- function(x,
                                        y,
                                        col = .crs_plot_color("data_overlay"),
                                        pch = .crs_plot_pch("data_overlay"),
                                        cex = .crs_plot_cex("data_overlay"),
                                        ...) {
  if (is.null(x) || is.null(y) || is.factor(x) || is.ordered(x))
    return(invisible(FALSE))
  ok <- is.finite(x) & is.finite(y)
  if (!any(ok)) return(invisible(FALSE))
  graphics::points(x[ok], y[ok], col = col, pch = pch, cex = cex, ...)
  invisible(TRUE)
}

.crs_plot_overlay_points_factor <- function(x,
                                            y,
                                            col = .crs_plot_color("data_overlay"),
                                            pch = .crs_plot_pch("data_overlay"),
                                            cex = .crs_plot_cex("data_overlay"),
                                            ...) {
  if (is.null(x) || is.null(y) || !(is.factor(x) || is.ordered(x)))
    return(invisible(FALSE))
  ok <- !is.na(x) & is.finite(y)
  if (!any(ok)) return(invisible(FALSE))
  graphics::points(x[ok], y[ok], col = col, pch = pch, cex = cex, ...)
  invisible(TRUE)
}

.crs_plot_draw_rug_1d <- function(x, col = .crs_plot_color("support"),
                                  quiet = TRUE, ...) {
  if (is.null(x) || is.factor(x) || is.ordered(x))
    return(invisible(FALSE))
  x <- as.vector(x)
  ok <- is.finite(x)
  if (!any(ok)) return(invisible(FALSE))
  graphics::rug(x[ok], col = col, quiet = quiet, ...)
  invisible(TRUE)
}

.crs_plot_overlay_points_persp <- function(x1,
                                           x2,
                                           y,
                                           persp.mat,
                                           col = .crs_plot_color("data_overlay"),
                                           pch = .crs_plot_pch("data_overlay"),
                                           cex = .crs_plot_cex("data_overlay"),
                                           ...) {
  if (is.null(x1) || is.null(x2) || is.null(y) || is.null(persp.mat))
    return(invisible(FALSE))
  ok <- is.finite(x1) & is.finite(x2) & is.finite(y)
  if (!any(ok)) return(invisible(FALSE))
  xy <- grDevices::trans3d(x1[ok], x2[ok], y[ok], persp.mat)
  graphics::points(xy$x, xy$y, col = col, pch = pch, cex = cex, ...)
  invisible(TRUE)
}

.crs_plot_draw_floor_rug_persp <- function(x1,
                                           x2,
                                           zlim,
                                           persp.mat,
                                           col = .crs_plot_color("support_floor"),
                                           lwd = .crs_plot_lwd("support"),
                                           ...) {
  if (is.null(x1) || is.null(x2) || is.null(zlim) || is.null(persp.mat))
    return(invisible(FALSE))
  ok <- is.finite(x1) & is.finite(x2)
  if (!any(ok)) return(invisible(FALSE))
  zlim <- as.double(zlim)
  if (length(zlim) < 2L || !all(is.finite(zlim))) return(invisible(FALSE))
  zfloor <- zlim[1L]
  zspan <- diff(range(zlim))
  if (!is.finite(zspan) || zspan <= 0) zspan <- max(1, abs(zfloor))
  ztop <- zfloor + 0.02 * zspan
  lower <- grDevices::trans3d(x1[ok], x2[ok],
                              rep.int(zfloor, sum(ok)), persp.mat)
  upper <- grDevices::trans3d(x1[ok], x2[ok],
                              rep.int(ztop, sum(ok)), persp.mat)
  graphics::segments(lower$x, lower$y, upper$x, upper$y,
                     col = col, lwd = lwd, ...)
  invisible(TRUE)
}

.crs_plot_draw_floor_rug_rgl <- function(x1,
                                         x2,
                                         zlim,
                                         color = .crs_plot_color("support_floor"),
                                         alpha = 0.40,
                                         lwd = .crs_plot_lwd("support"),
                                         ...) {
  if (is.null(x1) || is.null(x2) || is.null(zlim))
    return(invisible(FALSE))
  ok <- is.finite(x1) & is.finite(x2)
  if (!any(ok)) return(invisible(FALSE))
  zlim <- as.double(zlim)
  if (length(zlim) < 2L || !all(is.finite(zlim))) return(invisible(FALSE))
  zfloor <- zlim[1L]
  zspan <- diff(range(zlim))
  if (!is.finite(zspan) || zspan <= 0) zspan <- max(1, abs(zfloor))
  ztop <- zfloor + 0.02 * zspan
  rgl::segments3d(
    x = rbind(x1[ok], x1[ok]),
    y = rbind(x2[ok], x2[ok]),
    z = rbind(rep.int(zfloor, sum(ok)), rep.int(ztop, sum(ok))),
    color = color,
    alpha = alpha,
    lwd = lwd,
    ...
  )
  invisible(TRUE)
}

.crs_plot_draw_box_grid_persp <- function(xlim,
                                          ylim,
                                          zlim,
                                          persp.mat,
                                          col = .crs_plot_color("support_grid"),
                                          lwd = .crs_plot_lwd("support_grid")) {
  if (is.null(xlim) || is.null(ylim) || is.null(zlim) || is.null(persp.mat))
    return(invisible(FALSE))
  xlim <- as.double(xlim)
  ylim <- as.double(ylim)
  zlim <- as.double(zlim)
  if (length(xlim) < 2L || length(ylim) < 2L || length(zlim) < 2L)
    return(invisible(FALSE))
  if (!all(is.finite(c(xlim, ylim, zlim))))
    return(invisible(FALSE))

  x.at <- pretty(xlim, n = 5L)
  y.at <- pretty(ylim, n = 5L)
  z.at <- pretty(zlim, n = 5L)
  x.at <- x.at[x.at >= min(xlim) & x.at <= max(xlim)]
  y.at <- y.at[y.at >= min(ylim) & y.at <= max(ylim)]
  z.at <- z.at[z.at >= min(zlim) & z.at <= max(zlim)]

  draw_segment <- function(x0, y0, z0, x1, y1, z1) {
    pts <- grDevices::trans3d(c(x0, x1), c(y0, y1), c(z0, z1),
                              persp.mat)
    graphics::segments(pts$x[1L], pts$y[1L], pts$x[2L], pts$y[2L],
                       col = col, lwd = lwd)
  }

  xmin <- min(xlim); xmax <- max(xlim)
  ymin <- min(ylim); ymax <- max(ylim)
  zmin <- min(zlim); zmax <- max(zlim)
  for (yy in y.at) draw_segment(xmin, yy, zmin, xmax, yy, zmin)
  for (xx in x.at) draw_segment(xx, ymin, zmin, xx, ymax, zmin)
  for (yy in y.at) draw_segment(xmin, yy, zmin, xmin, yy, zmax)
  for (zz in z.at) draw_segment(xmin, ymin, zz, xmin, ymax, zz)
  for (xx in x.at) draw_segment(xx, ymax, zmin, xx, ymax, zmax)
  for (zz in z.at) draw_segment(xmin, ymax, zz, xmax, ymax, zz)
  invisible(TRUE)
}

.crs_plot_wireframe_templates_persp <- function(x, y) {
  x <- as.double(x)
  y <- as.double(y)
  nx <- length(x)
  ny <- length(y)
  list(
    row_x = matrix(rep(x, ny), nrow = nx, ncol = ny),
    row_y = matrix(rep(y, each = nx), nrow = nx, ncol = ny),
    col_x = matrix(rep(x, each = ny), nrow = ny, ncol = nx),
    col_y = matrix(rep(y, nx), nrow = ny, ncol = nx),
    nx = nx,
    ny = ny
  )
}

.crs_plot_project_wire_matrix_persp <- function(xmat, ymat, zmat, pmat) {
  tr <- cbind(as.vector(xmat), as.vector(ymat), as.vector(zmat), 1) %*% pmat
  list(
    x = matrix(tr[, 1L] / tr[, 4L], nrow = nrow(xmat), ncol = ncol(xmat)),
    y = matrix(tr[, 2L] / tr[, 4L], nrow = nrow(ymat), ncol = ncol(ymat))
  )
}

.crs_plot_draw_wire_surface_persp <- function(templates,
                                              zmat,
                                              persp.mat,
                                              col,
                                              lwd = 1) {
  if (is.null(zmat) || !any(is.finite(zmat)))
    return(invisible(FALSE))
  row_pts <- .crs_plot_project_wire_matrix_persp(
    xmat = templates$row_x,
    ymat = templates$row_y,
    zmat = zmat,
    pmat = persp.mat
  )
  for (j in seq_len(templates$ny))
    graphics::lines(row_pts$x[, j], row_pts$y[, j], col = col, lwd = lwd)

  col_pts <- .crs_plot_project_wire_matrix_persp(
    xmat = templates$col_x,
    ymat = templates$col_y,
    zmat = t(zmat),
    pmat = persp.mat
  )
  for (i in seq_len(templates$nx))
    graphics::lines(col_pts$x[, i], col_pts$y[, i], col = col, lwd = lwd)
  invisible(TRUE)
}

.crs_plot_all_band_colors <- function() {
  .crs_plot_viridis_at(
    c(pointwise = 0.18, simultaneous = 0.50, bonferroni = 0.78)
  )
}

.crs_plot_all_band_alpha <- function() {
  c(pointwise = 0.14,
    simultaneous = 0.10,
    bonferroni = 0.08)
}

.crs_plot_all_band_legend <- function(legend = TRUE, where = "topright",
                                      lty = .crs_plot_lty("solid"),
                                      lwd = .crs_plot_lwd("band_all_1d"),
                                      ...) {
  if (is.null(legend)) legend <- TRUE
  if (identical(legend, FALSE)) return(invisible(FALSE))

  cols <- .crs_plot_all_band_colors()
  args <- list(
    x = where,
    legend = c("Pointwise", "Simultaneous", "Bonferroni"),
    col = unname(cols[c("pointwise", "simultaneous", "bonferroni")]),
    lty = lty,
    lwd = lwd,
    bty = "n",
    cex = .crs_plot_cex("legend")
  )

  if (is.character(legend) && length(legend) == 1L) {
    args$x <- legend
  } else if (is.list(legend)) {
    for (nm in names(legend)) args[[nm]] <- legend[[nm]]
  } else if (!is.logical(legend) || length(legend) != 1L || is.na(legend)) {
    stop("legend must be TRUE, FALSE, a legend position, or a list",
         call. = FALSE)
  }

  do.call(graphics::legend, args)
  invisible(TRUE)
}

.crs_plot_merge_rgl_legend_control <- function(legend3d.args,
                                               legend = TRUE) {
  legend.value <- if (is.list(legend) &&
                      any(names(legend) %in% c("tau", "bands"))) {
    if (!is.null(legend$bands)) legend$bands else TRUE
  } else {
    legend
  }

  if (is.null(legend.value))
    return(.crs_plot_merge_override_args(legend3d.args, list(plot = FALSE)))

  if (is.logical(legend.value)) {
    if (length(legend.value) != 1L)
      stop("legend must be TRUE/FALSE, NULL, NA, a legend position string, or a list of graphics::legend arguments",
           call. = FALSE)
    if (is.na(legend.value) || !isTRUE(legend.value))
      return(.crs_plot_merge_override_args(legend3d.args,
                                           list(plot = FALSE)))
    return(legend3d.args)
  }

  if (is.character(legend.value) && length(legend.value) == 1L &&
      !is.na(legend.value))
    return(.crs_plot_merge_override_args(list(x = legend.value),
                                         legend3d.args))

  if (is.list(legend.value)) {
    show <- legend.value$show
    if (!is.null(show)) {
      if (!is.logical(show) || length(show) != 1L)
        stop("legend$show must be TRUE or FALSE", call. = FALSE)
      legend.value$show <- NULL
      if (is.na(show) || !isTRUE(show))
        return(.crs_plot_merge_override_args(legend3d.args,
                                             list(plot = FALSE)))
    }
    return(.crs_plot_merge_override_args(legend.value, legend3d.args))
  }

  stop("legend must be TRUE/FALSE, NULL, NA, a legend position string, or a list of graphics::legend arguments",
       call. = FALSE)
}

.crs_plot_draw_error_wireframes_persp <- function(x,
                                                  y,
                                                  persp.mat,
                                                  plot.errors.type,
                                                  lerr = NULL,
                                                  herr = NULL,
                                                  lerr.all = NULL,
                                                  herr.all = NULL,
                                                  lwd = graphics::par()$lwd) {
  templates <- .crs_plot_wireframe_templates_persp(x = x, y = y)
  if (identical(plot.errors.type, "all") &&
      !is.null(lerr.all) && !is.null(herr.all)) {
    band.cols <- .crs_plot_all_band_colors()
    for (bn in c("pointwise", "simultaneous", "bonferroni")) {
      col <- grDevices::adjustcolor(band.cols[[bn]], alpha.f = 0.50)
      .crs_plot_draw_wire_surface_persp(
        templates, lerr.all[[bn]], persp.mat, col = col,
        lwd = .crs_plot_lwd("band_all_surface", lwd)
      )
      .crs_plot_draw_wire_surface_persp(
        templates, herr.all[[bn]], persp.mat, col = col,
        lwd = .crs_plot_lwd("band_all_surface", lwd)
      )
    }
    return(invisible(TRUE))
  }
  wire.col <- grDevices::adjustcolor(.crs_plot_color("primary"),
                                     alpha.f = 0.45)
  .crs_plot_draw_wire_surface_persp(
    templates, lerr, persp.mat, col = wire.col,
    lwd = .crs_plot_lwd("interval_surface", lwd)
  )
  .crs_plot_draw_wire_surface_persp(
    templates, herr, persp.mat, col = wire.col,
    lwd = .crs_plot_lwd("interval_surface", lwd)
  )
  invisible(TRUE)
}

.crs_plot_error_surfaces_rgl <- function(x,
                                         y,
                                         plot.errors.type,
                                         lerr = NULL,
                                         herr = NULL,
                                         lerr.all = NULL,
                                         herr.all = NULL,
                                         surface3d.args = list(),
                                         legend3d.args = list(),
                                         ...) {
  draw_one <- function(z, color) {
    if (is.null(z) || !any(is.finite(z)))
      return(invisible(FALSE))
    surf.args <- .crs_plot_merge_override_args(
      list(x = x, y = y, z = z, color = color, alpha = 0.20,
           front = "lines", back = "lines", lit = FALSE),
      .crs_plot_merge_override_args(surface3d.args, list(...))
    )
    do.call(rgl::surface3d, surf.args)
    invisible(TRUE)
  }
  if (identical(plot.errors.type, "all") &&
      !is.null(lerr.all) && !is.null(herr.all)) {
    band.cols <- .crs_plot_all_band_colors()
    band.alpha <- .crs_plot_all_band_alpha()
    drawn.bands <- character(0L)
    for (bn in c("pointwise", "simultaneous", "bonferroni")) {
      col <- grDevices::adjustcolor(band.cols[[bn]],
                                    alpha.f = band.alpha[[bn]])
      drawn.lower <- draw_one(lerr.all[[bn]], col)
      drawn.upper <- draw_one(herr.all[[bn]], col)
      if (isTRUE(drawn.lower) || isTRUE(drawn.upper))
        drawn.bands <- c(drawn.bands, bn)
    }
    if (!length(drawn.bands))
      return(invisible(FALSE))
    legend3d.call <- .crs_plot_merge_override_args(
      list(
        "topright",
        legend = c(pointwise = "Pointwise",
                   simultaneous = "Simultaneous",
                   bonferroni = "Bonferroni")[drawn.bands],
        col = unname(band.cols[drawn.bands]),
        lty = .crs_plot_lty("solid"),
        lwd = .crs_plot_lwd("band_all_surface"),
        cex = .crs_plot_cex("legend"),
        bg = .crs_plot_color("legend_bg"),
        bty = "n"
      ),
      legend3d.args
    )
    do.call(rgl::legend3d, legend3d.call)
    return(invisible(TRUE))
  }
  col <- grDevices::adjustcolor(.crs_plot_color("primary"), alpha.f = 0.45)
  draw_one(lerr, col)
  draw_one(herr, col)
  invisible(TRUE)
}

.crs_plot_stop_unused_args <- function(bad, allowed) {
  bad <- unique(bad[nzchar(bad)])
  if (!length(bad)) return(invisible(NULL))
  msg <- if (length(bad) == 1L) {
    sprintf("unused plot argument: %s", bad)
  } else {
    sprintf("unused plot arguments: %s", paste(bad, collapse = ", "))
  }
  close <- vapply(bad, function(x) {
    d <- utils::adist(x, allowed, partial = FALSE, ignore.case = FALSE)
    if (!length(d)) return(NA_character_)
    i <- which.min(d)
    if (is.finite(d[i]) && d[i] <= max(2L, floor(nchar(x) / 3L)))
      allowed[i]
    else
      NA_character_
  }, character(1L))
  close <- unique(stats::na.omit(close))
  if (length(close))
    msg <- paste0(msg, "; did you mean ",
                  paste(close, collapse = " or "), "?")
  stop(msg, call. = FALSE)
}

.crs_plot_graphics_arg_names <- function() {
  unique(c(
    setdiff(names(formals(graphics::plot.default)), c("x", "y", "...")),
    names(graphics::par(no.readonly = TRUE)),
    "panel.first", "panel.last", "zlab", "zlim", "theta", "phi", "border",
    "view", "type", "lty", "lwd", "col", "pch", "cex", "main", "sub",
    "xlab", "ylab", "xlim", "ylim"
  ))
}

.crs_plot_canonical_arg_names <- function() {
  c("errors", "band", "alpha", "bootstrap", "B", "center",
    "output", "data_overlay", "data_rug", "layout", "legend",
    "factor_boxplot", "boxplot_outliers",
    "gradient", "gradients", "gradient.order", "gradient_order",
    "common_scale",
    "renderer", "neval", "perspective", "view", "behavior",
    "boot_control", "grid_control", "render_control")
}

.crs_plot_legacy_arg_names <- function() {
  c("plot.errors.method", "plot.errors.type", "plot.errors.alpha",
    "plot.errors.boot.method", "plot.errors.boot.num",
    "plot.errors.boot.nonfixed", "plot.errors.boot.wild",
    "plot.errors.boot.blocklen", "plot.errors.center",
    "plot.errors.style", "plot.errors.bar", "plot.errors.bar.num",
    "plot.behavior", "plot.data.overlay", "plot.rug", "plot.par.mfrow",
    "plot.bxp", "plot.bxp.out", "num.eval", "persp", "xtrim", "xq",
    "common.scale", "display.nomad.progress", "display.warnings")
}

.crs_plot_validate_public_dots <- function(dots_call,
                                           context = "plot",
                                           extra = character()) {
  nms <- .crs_plot_dot_names(dots_call)
  if (any(!nzchar(nms)))
    stop(sprintf("unnamed plot arguments are not supported for %s", context),
         call. = FALSE)
  if ("intervals" %in% nms)
    stop("unused plot argument: intervals; did you mean errors?",
         call. = FALSE)
  if ("boot" %in% nms)
    stop("unused plot argument: boot; did you mean bootstrap?",
         call. = FALSE)
  if ("bands" %in% nms)
    stop("unused plot argument: bands; did you mean band?",
         call. = FALSE)
  if ("ci" %in% nms)
    stop("unused plot argument: ci; use errors instead",
         call. = FALSE)
  if ("deriv" %in% nms)
    stop("unused plot argument: deriv; use gradients and gradient_order instead",
         call. = FALSE)
  if ("mean" %in% nms)
    stop("unused plot argument: mean; plot.crs displays fitted functions by default",
         call. = FALSE)
  if ("plot.view" %in% nms)
    stop("unused plot argument: plot.view; CRS plot methods now use the NP plot interface",
         call. = FALSE)
  legacy <- intersect(nms, .crs_plot_legacy_arg_names())
  if (length(legacy)) {
    stop(sprintf("unused plot argument: %s; use NP-style plot arguments",
                 legacy[1L]),
         call. = FALSE)
  }
  allowed <- unique(c(.crs_plot_graphics_arg_names(),
                      .crs_plot_canonical_arg_names(),
                      extra))
  bad <- nms[nzchar(nms) & !(nms %in% allowed) &
               !grepl("^(plot|lines|points|persp|rgl\\.[A-Za-z0-9_]+)\\.",
                      nms)]
  .crs_plot_stop_unused_args(bad, allowed)
  invisible(TRUE)
}

.crs_plot_reject_curve_controls <- function(dots,
                                            context = "plot",
                                            allow.data.overlay = FALSE) {
  if (is.null(dots) || !length(dots)) return(invisible(TRUE))
  nms <- names(dots)
  if (is.null(nms)) return(invisible(TRUE))

  unsupported <- c("errors", "band", "alpha", "bootstrap", "B", "center",
                   "data_rug", "layout", "legend",
                   "factor_boxplot", "boxplot_outliers",
                   "gradient", "gradients", "gradient.order",
                   "gradient_order", "common_scale",
                   "renderer", "neval", "perspective", "view",
                   "boot_control", "grid_control", "render_control")
  if (!isTRUE(allow.data.overlay))
    unsupported <- c(unsupported, "data_overlay")
  bad <- intersect(nms, unsupported)
  if (length(bad)) {
    stop(sprintf("%s does not support plot argument %s; this curve route supports output/behavior plus ordinary graphics arguments",
                 context, bad[1L]),
         call. = FALSE)
  }
  invisible(TRUE)
}

.crs_plot_response_label <- function(object, fallback = "Conditional Mean") {
  if (!is.null(object$formula)) {
    vars <- all.vars(object$formula)
    if (length(vars) >= 1L && nzchar(vars[1L]))
      return(vars[1L])
  }
  if (!is.null(object$terms)) {
    vars <- attr(object$terms, "variables")
    if (length(vars) >= 2L) {
      lab <- deparse(vars[[2L]], width.cutoff = 500L)
      if (length(lab) && nzchar(lab[1L]))
        return(lab[1L])
    }
  }
  fallback
}

.crs_plot_draw_factor_fit <- function(x, y, col = graphics::par()$col,
                                      lty = .crs_plot_lty("interval"),
                                      lwd = graphics::par()$lwd,
                                      pch = .crs_plot_pch("fit"),
                                      cex = 1, bg = NULL) {
  l.f <- rep(x, each = 3L)
  l.f[3L * seq_along(x)] <- NA
  l.y <- unlist(lapply(y, function(p) c(0, p, NA)), use.names = FALSE)
  graphics::lines(x = l.f, y = l.y, col = col, lty = lty, lwd = lwd)
  point.args <- list(x = x, y = y, col = col, pch = pch, cex = cex)
  if (!is.null(bg)) point.args$bg <- bg
  do.call(graphics::points, point.args)
  invisible(TRUE)
}

.crs_plot_draw_interval_bars <- function(x, lower, upper, col,
                                         lty = .crs_plot_lty("interval"),
                                         lwd = graphics::par()$lwd,
                                         cap = 0.08) {
  ok <- is.finite(as.numeric(x)) & is.finite(lower) & is.finite(upper)
  if (!any(ok)) return(invisible(FALSE))
  xx <- as.numeric(x[ok])
  lower <- lower[ok]
  upper <- upper[ok]
  graphics::segments(xx, lower, xx, upper, col = col, lty = lty, lwd = lwd)
  graphics::segments(xx - cap, lower, xx + cap, lower,
                     col = col, lty = lty, lwd = lwd)
  graphics::segments(xx - cap, upper, xx + cap, upper,
                     col = col, lty = lty, lwd = lwd)
  invisible(TRUE)
}

.crs_plot_set_normalized_arg <- function(dots, public, internal, value) {
  same <- !is.null(dots[[internal]]) &&
    isTRUE(all.equal(dots[[internal]], value, check.attributes = FALSE))
  if (!is.null(dots[[internal]]) && !same)
    stop(sprintf("conflicting plot arguments: %s and %s specify different values",
                 public, internal),
         call. = FALSE)
  dots[[internal]] <- value
  dots
}

.crs_plot_match_layout <- function(value) {
  if (is.logical(value)) {
    if (length(value) != 1L || is.na(value))
      stop("layout must be TRUE/FALSE or one of \"auto\", \"current\"",
           call. = FALSE)
    return(isTRUE(value))
  }
  layout <- .crs_plot_scalar_match(value, c("auto", "current"), "layout")
  identical(layout, "auto")
}

.crs_plot_normalize_public_dots <- function(dots, context = "plot") {
  if (is.null(dots)) dots <- list()
  supplied <- names(dots)
  has <- function(nm) !is.null(dots[[nm]])

  if (has("intervals"))
    stop("unused plot argument: intervals; did you mean errors?",
         call. = FALSE)
  if (has("boot"))
    stop("unused plot argument: boot; did you mean bootstrap?",
         call. = FALSE)
  if (has("bands"))
    stop("unused plot argument: bands; did you mean band?",
         call. = FALSE)

  if (has("errors")) {
    errors <- .crs_plot_scalar_match(dots$errors,
                                     c("none", "bootstrap", "asymptotic"),
                                     "errors")
    dots$errors <- NULL
    dots <- .crs_plot_set_normalized_arg(dots, "errors",
                                         "plot.errors.method", errors)
  }
  if (has("band")) {
    band <- .crs_plot_scalar_match(dots$band,
                                   c("pmzsd", "pointwise", "bonferroni",
                                     "simultaneous", "all"),
                                   "band")
    dots$band <- NULL
    if (identical(band, "pmzsd")) band <- "standard"
    dots <- .crs_plot_set_normalized_arg(dots, "band",
                                         "plot.errors.type", band)
  }
  if (has("alpha")) {
    alpha <- dots$alpha
    if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
        alpha <= 0 || alpha >= 0.5)
      stop("alpha must lie in (0, 0.5)", call. = FALSE)
    dots$alpha <- NULL
    dots <- .crs_plot_set_normalized_arg(dots, "alpha",
                                         "plot.errors.alpha", alpha)
  }
  if (has("bootstrap")) {
    if (is.list(dots$bootstrap))
      stop("unused plot argument: bootstrap list; use scalar bootstrap, B, and boot_control",
           call. = FALSE)
    bootstrap <- .crs_plot_scalar_match(dots$bootstrap,
                                        c("wild", "inid", "fixed", "geom"),
                                        "bootstrap")
    dots$bootstrap <- NULL
    dots <- .crs_plot_set_normalized_arg(dots, "bootstrap",
                                         "plot.errors.boot.method", bootstrap)
  }
  if (has("B")) {
    B <- dots$B
    if (!is.numeric(B) || length(B) != 1L || is.na(B) || B < 1)
      stop("B must be a positive numeric scalar", call. = FALSE)
    dots$B <- NULL
    dots <- .crs_plot_set_normalized_arg(dots, "B",
                                         "plot.errors.boot.num",
                                         as.integer(B))
  }
  if (has("center")) {
    center <- .crs_plot_scalar_match(dots$center,
                                     c("estimate", "bias-corrected"),
                                     "center")
    dots$center <- NULL
    dots <- .crs_plot_set_normalized_arg(dots, "center",
                                         "plot.errors.center", center)
  }
  if (has("output")) {
    output <- .crs_plot_scalar_match(dots$output,
                                     c("plot", "data", "both", "plot-data"),
                                     "output")
    if (identical(output, "both")) output <- "plot-data"
    dots$output <- NULL
    dots <- .crs_plot_set_normalized_arg(dots, "output",
                                         "plot.behavior", output)
  }
  if (has("behavior")) {
    behavior <- .crs_plot_scalar_match(dots$behavior,
                                       c("plot", "plot-data", "data"),
                                       "behavior")
    dots$behavior <- NULL
    dots <- .crs_plot_set_normalized_arg(dots, "behavior",
                                         "plot.behavior", behavior)
  }
  if (has("data_overlay")) {
    data_overlay <- .crs_plot_match_flag(dots$data_overlay, "data_overlay")
    dots$data_overlay <- NULL
    dots <- .crs_plot_set_normalized_arg(dots, "data_overlay",
                                         "plot.data.overlay", data_overlay)
  }
  if (has("data_rug")) {
    data_rug <- .crs_plot_match_flag(dots$data_rug, "data_rug")
    dots$data_rug <- NULL
    dots <- .crs_plot_set_normalized_arg(dots, "data_rug",
                                         "plot.rug", data_rug)
  }
  if (has("layout")) {
    layout <- .crs_plot_match_layout(dots$layout)
    dots$layout <- NULL
    dots <- .crs_plot_set_normalized_arg(dots, "layout",
                                         "plot.par.mfrow", layout)
  }
  if (has("factor_boxplot")) {
    factor_boxplot <- .crs_plot_match_flag(dots$factor_boxplot,
                                           "factor_boxplot")
    dots$factor_boxplot <- NULL
    dots <- .crs_plot_set_normalized_arg(dots, "factor_boxplot",
                                         "plot.bxp", factor_boxplot)
  }
  if (has("boxplot_outliers")) {
    boxplot_outliers <- .crs_plot_match_flag(dots$boxplot_outliers,
                                             "boxplot_outliers")
    dots$boxplot_outliers <- NULL
    dots <- .crs_plot_set_normalized_arg(dots, "boxplot_outliers",
                                         "plot.bxp.out", boxplot_outliers)
  }
  if (has("neval")) {
    neval <- dots$neval
    if (!is.numeric(neval) || length(neval) != 1L || is.na(neval) ||
        neval < 2)
      stop("neval must be a numeric scalar >= 2", call. = FALSE)
    dots$neval <- NULL
    dots <- .crs_plot_set_normalized_arg(dots, "neval", "num.eval",
                                         as.integer(neval))
  }
  if (has("gradient_order")) {
    gradient_order <- dots$gradient_order
    if (!is.numeric(gradient_order) || any(is.na(gradient_order)) ||
        any(gradient_order < 1L))
      stop("gradient_order must contain positive numeric values",
           call. = FALSE)
    dots$gradient_order <- NULL
    dots <- .crs_plot_set_normalized_arg(dots, "gradient_order",
                                         "gradient.order",
                                         as.integer(gradient_order))
  }
  if (has("gradient")) {
    gradient <- .crs_plot_match_flag(dots$gradient, "gradient")
    dots$gradient <- NULL
    dots <- .crs_plot_set_normalized_arg(dots, "gradient",
                                         "gradients", gradient)
  }
  if (has("common_scale")) {
    common_scale <- .crs_plot_match_flag(dots$common_scale, "common_scale")
    dots$common_scale <- NULL
    dots <- .crs_plot_set_normalized_arg(dots, "common_scale",
                                         "common.scale", common_scale)
  }
  if (has("perspective")) {
    perspective <- .crs_plot_match_flag(dots$perspective, "perspective")
    dots$perspective <- NULL
    dots <- .crs_plot_set_normalized_arg(dots, "perspective",
                                         "perspective", perspective)
  }
  if (has("renderer")) {
    renderer <- .crs_plot_scalar_match(dots$renderer, c("base", "rgl"),
                                       "renderer")
    dots$renderer <- NULL
    dots <- .crs_plot_set_normalized_arg(dots, "renderer", "renderer",
                                         renderer)
  }
  if (has("boot_control")) {
    if (!inherits(dots$boot_control, "np_boot_control") &&
        !inherits(dots$boot_control, "crs_boot_control"))
      stop("boot_control must be created by np_boot_control()",
           call. = FALSE)
    ctrl <- dots$boot_control
    dots$boot_control <- NULL
    dots <- .crs_plot_set_normalized_arg(dots, "boot_control$nonfixed",
                                         "plot.errors.boot.nonfixed",
                                         ctrl$nonfixed)
    dots <- .crs_plot_set_normalized_arg(dots, "boot_control$wild",
                                         "plot.errors.boot.wild", ctrl$wild)
    if (!is.null(ctrl$blocklen))
      dots <- .crs_plot_set_normalized_arg(dots, "boot_control$blocklen",
                                           "plot.errors.boot.blocklen",
                                           ctrl$blocklen)
  }
  if (has("grid_control")) {
    if (!inherits(dots$grid_control, "np_grid_control") &&
        !inherits(dots$grid_control, "crs_grid_control"))
      stop("grid_control must be created by np_grid_control()",
           call. = FALSE)
    ctrl <- dots$grid_control
    dots$grid_control <- NULL
    if (!is.null(ctrl$xtrim))
      dots <- .crs_plot_set_normalized_arg(dots, "grid_control$xtrim",
                                           "xtrim", ctrl$xtrim)
    if (!is.null(ctrl$xq))
      dots <- .crs_plot_set_normalized_arg(dots, "grid_control$xq",
                                           "xq", ctrl$xq)
    if (!is.null(ctrl$slices))
      stop(sprintf("grid_control$slices is not yet supported for %s", context),
           call. = FALSE)
  }
  if (has("render_control")) {
    if (!inherits(dots$render_control, "np_render_control") &&
        !inherits(dots$render_control, "crs_render_control"))
      stop("render_control must be created by np_render_control()",
           call. = FALSE)
    ctrl <- dots$render_control
    dots$render_control <- NULL
    dots <- .crs_plot_set_normalized_arg(dots, "render_control$style",
                                         "plot.errors.style", ctrl$style)
    dots <- .crs_plot_set_normalized_arg(dots, "render_control$bar",
                                         "plot.errors.bar", ctrl$bar)
    if (!is.null(ctrl$bar_num))
      dots <- .crs_plot_set_normalized_arg(dots, "render_control$bar_num",
                                           "plot.errors.bar.num",
                                           ctrl$bar_num)
  }

  method <- if (!is.null(dots$plot.errors.method)) {
    as.character(dots$plot.errors.method)[1L]
  } else {
    "none"
  }
  boot.only <- c("plot.errors.boot.method", "plot.errors.boot.num",
                 "plot.errors.boot.nonfixed", "plot.errors.boot.wild",
                 "plot.errors.boot.blocklen")
  boot.supplied <- any(c("bootstrap", "B", "boot_control", "center",
                         "plot.errors.center", boot.only) %in% supplied)
  if (!identical(method, "bootstrap") && boot.supplied)
    stop("bootstrap controls require errors = \"bootstrap\"",
         call. = FALSE)

  dots
}
