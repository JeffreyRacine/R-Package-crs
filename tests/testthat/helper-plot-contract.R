crs_plot_contract_cases <- function() {
  list(
    canonical = c("errors", "band", "alpha", "bootstrap", "B", "center",
                  "output", "data_overlay", "data_rug", "layout", "legend",
                  "factor_boxplot", "boxplot_outliers", "gradient",
                  "gradients", "gradient_order", "common_scale", "xtrim",
                  "xq", "renderer", "neval", "perspective", "view", "behavior",
                  "boot_control", "grid_control", "render_control"),
    stale = c("ci", "deriv", "mean", "plot.view", "intervals", "boot",
              "bands", "plot.errors.method", "plot.errors.type",
              "plot.errors.alpha", "plot.errors.boot.method",
              "plot.errors.boot.num", "plot.errors.boot.nonfixed",
              "plot.errors.boot.wild", "plot.errors.boot.blocklen",
              "plot.errors.center", "plot.errors.style", "plot.errors.bar",
              "plot.errors.bar.num", "plot.behavior", "plot.data.overlay",
              "plot.rug", "plot.par.mfrow", "plot.bxp", "plot.bxp.out",
              "num.eval", "persp", "common.scale",
              "display.nomad.progress", "display.warnings")
  )
}

crs_plot_pdf <- function(expr) {
  pdf.file <- tempfile(fileext = ".pdf")
  grDevices::pdf(pdf.file)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
}
