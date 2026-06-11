## Internal CRS plot helpers copied/adapted from the modern np plot layer.
##
## This file intentionally starts small. It provides the normalization and
## rendering-adjacent primitives needed by shadow plot engines without changing
## any public plot method defaults.

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

.crs_boot_control <- function(nonfixed = c("exact", "frozen"),
                              wild = c("rademacher", "mammen"),
                              blocklen = NULL) {
  nonfixed <- match.arg(nonfixed)
  wild <- match.arg(wild)
  if (!is.null(blocklen) &&
      (!is.numeric(blocklen) || length(blocklen) != 1L ||
       is.na(blocklen) || blocklen <= 0))
    stop("blocklen must be a positive numeric scalar", call. = FALSE)
  structure(list(nonfixed = nonfixed, wild = wild, blocklen = blocklen),
            class = "crs_boot_control")
}

.crs_grid_control <- function(xtrim = NULL, xq = NULL, slices = NULL) {
  if (!is.null(xtrim) &&
      (!is.numeric(xtrim) || length(xtrim) != 2L ||
       any(is.na(xtrim)) || any(xtrim < 0) || any(xtrim > 1) ||
       xtrim[1L] >= xtrim[2L]))
    stop("xtrim must be a numeric length-two vector with 0 <= xtrim[1] < xtrim[2] <= 1",
         call. = FALSE)
  structure(list(xtrim = xtrim, xq = xq, slices = slices),
            class = "crs_grid_control")
}

.crs_render_control <- function(style = c("band", "bar"),
                                bar = c("|", "I"),
                                bar_num = NULL) {
  style <- match.arg(style)
  bar <- match.arg(bar)
  if (!is.null(bar_num) &&
      (!is.numeric(bar_num) || length(bar_num) != 1L ||
       is.na(bar_num) || bar_num < 1))
    stop("bar_num must be a positive numeric scalar", call. = FALSE)
  structure(list(style = style, bar = bar, bar_num = bar_num),
            class = "crs_render_control")
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

.crs_plot_color <- function(role, alpha = NULL) {
  roles <- list(
    primary = "#0072B2",
    secondary = "#D55E00",
    fit = "#0072B2",
    interval = "#D55E00",
    interval_context = "#D55E00",
    data_overlay = "#000000",
    support = "#999999",
    surface_border = "#3A3A3A",
    context_wire = "#D55E00",
    context_border = "#555555",
    legend_bg = "#FFFFFF"
  )
  col <- roles[[role]]
  if (is.null(col)) col <- role
  if (!is.null(alpha)) col <- grDevices::adjustcolor(col, alpha.f = alpha)
  col
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

.crs_plot_user_args <- function(dots, type) {
  if (is.null(dots) || !length(dots)) return(list())
  prefix <- paste0(type, ".")
  nms <- names(dots)
  if (is.null(nms)) return(list())

  direct <- dots[[type]]
  if (!is.null(direct)) {
    if (!is.list(direct))
      stop(sprintf("%s must be a list when supplied as a grouped plot argument",
                   type),
           call. = FALSE)
  } else {
    direct <- list()
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
    "gradient_order", "common_scale", "renderer", "neval", "perspective",
    "boot_control", "grid_control", "render_control")
}

.crs_plot_legacy_arg_names <- function() {
  c("plot.errors.method", "plot.errors.type", "plot.errors.alpha",
    "plot.errors.boot.method", "plot.errors.boot.num",
    "plot.errors.boot.nonfixed", "plot.errors.boot.wild",
    "plot.errors.boot.blocklen", "plot.errors.center",
    "plot.errors.style", "plot.errors.bar", "plot.errors.bar.num",
    "plot.behavior", "plot.data.overlay", "plot.rug", "plot.par.mfrow",
    "plot.bxp", "plot.bxp.out", "num.eval", "persp", "persp.rgl",
    "mean", "deriv", "ci", "xtrim", "xq", "common.scale", "plot.view",
    "display.nomad.progress", "display.warnings")
}

.crs_plot_validate_public_dots <- function(dots_call,
                                           context = "plot",
                                           extra = character()) {
  nms <- if (is.null(dots_call) || length(dots_call) == 0L) {
    character()
  } else {
    names(dots_call)
  }
  if (is.null(nms)) nms <- rep.int("", length(dots_call))
  allowed <- unique(c(.crs_plot_graphics_arg_names(),
                      .crs_plot_canonical_arg_names(),
                      .crs_plot_legacy_arg_names(),
                      extra))
  bad <- nms[nzchar(nms) & !(nms %in% allowed) &
               !grepl("^(plot|lines|points|persp|rgl\\.[A-Za-z0-9_]+)\\.",
                      nms)]
  .crs_plot_stop_unused_args(bad, allowed)
  invisible(TRUE)
}

.crs_plot_set_normalized_arg <- function(dots, public, internal, value) {
  if (!is.null(dots[[internal]]))
    stop(sprintf("cannot supply both %s and %s", public, internal),
         call. = FALSE)
  dots[[internal]] <- value
  dots
}

.crs_plot_normalize_public_dots <- function(dots, context = "plot") {
  if (is.null(dots)) dots <- list()
  supplied <- names(dots)
  has <- function(nm) !is.null(dots[[nm]])

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
                                   c("standard", "pointwise", "bonferroni",
                                     "simultaneous", "all"),
                                   "band")
    dots$band <- NULL
    dots <- .crs_plot_set_normalized_arg(dots, "band",
                                         "plot.errors.type", band)
  }
  if (has("alpha")) {
    alpha <- dots$alpha
    if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
        alpha <= 0 || alpha >= 1)
      stop("alpha must be a numeric scalar in (0, 1)", call. = FALSE)
    dots$alpha <- NULL
    dots <- .crs_plot_set_normalized_arg(dots, "alpha",
                                         "plot.errors.alpha", alpha)
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
  if (has("output")) {
    output <- .crs_plot_scalar_match(dots$output, c("plot", "plot-data", "data"),
                                     "output")
    dots$output <- NULL
    dots <- .crs_plot_set_normalized_arg(dots, "output",
                                         "plot.behavior", output)
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
  if (has("neval")) {
    neval <- dots$neval
    if (!is.numeric(neval) || length(neval) != 1L || is.na(neval) ||
        neval < 2)
      stop("neval must be a numeric scalar >= 2", call. = FALSE)
    dots$neval <- NULL
    dots <- .crs_plot_set_normalized_arg(dots, "neval", "num.eval",
                                         as.integer(neval))
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
    if (!inherits(dots$boot_control, "crs_boot_control"))
      stop("boot_control must be created by .crs_boot_control()",
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
    if (!inherits(dots$grid_control, "crs_grid_control"))
      stop("grid_control must be created by .crs_grid_control()",
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
    if (!inherits(dots$render_control, "crs_render_control"))
      stop("render_control must be created by .crs_render_control()",
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
