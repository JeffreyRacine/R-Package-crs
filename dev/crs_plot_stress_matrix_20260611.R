#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  i <- match(name, args)
  if (is.na(i) || i == length(args)) return(default)
  args[[i + 1L]]
}

repo <- normalizePath(get_arg("--repo", getwd()), mustWork = FALSE)
lib <- normalizePath(
  get_arg("--lib", "/Users/jracine/Development/tmp/crs_plot_contract_20260611/Rlib"),
  mustWork = FALSE
)
out_dir <- normalizePath(
  get_arg("--out-dir", tempfile("crs_plot_stress_")),
  mustWork = FALSE
)
B <- as.integer(get_arg("--B", "5"))
neval <- as.integer(get_arg("--neval", "7"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "warnings"), recursive = TRUE, showWarnings = FALSE)

if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(crs))

strip_time <- function(x) {
  if (is.list(x)) {
    x$timing <- NULL
    x$plot.timing <- NULL
    x$total.time <- NULL
    x <- lapply(x, strip_time)
  }
  x
}

has_bad_numeric <- function(x) {
  if (is.numeric(x)) return(any(is.nan(x) | is.infinite(x)))
  if (is.data.frame(x)) return(any(vapply(x, has_bad_numeric, logical(1L))))
  if (is.list(x)) return(any(vapply(strip_time(x), has_bad_numeric, logical(1L))))
  FALSE
}

run_pdf <- function(expr) {
  pdf_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(pdf_file)
  on.exit({
    grDevices::dev.off()
    unlink(pdf_file)
  }, add = TRUE)
  force(expr)
}

run_case <- function(name,
                     expected = c("ok", "error"),
                     pattern = "",
                     max_warnings = Inf,
                     expr) {
  expected <- match.arg(expected)
  warn_log <- character()
  err <- ""
  status <- "ok"
  value <- NULL
  t0 <- proc.time()[["elapsed"]]
  value <- tryCatch(
    withCallingHandlers(
      force(expr),
      warning = function(w) {
        warn_log <<- c(warn_log, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      status <<- "error"
      err <<- conditionMessage(e)
      NULL
    }
  )
  elapsed <- proc.time()[["elapsed"]] - t0
  nonfinite <- if (!is.null(value)) isTRUE(has_bad_numeric(value)) else FALSE
  warn_ok <- length(warn_log) <= max_warnings
  pass <- if (identical(expected, "ok")) {
    identical(status, "ok") && !isTRUE(nonfinite) && isTRUE(warn_ok)
  } else {
    identical(status, "error") &&
      (!nzchar(pattern) || grepl(pattern, err, fixed = FALSE))
  }
  if (length(warn_log)) {
    writeLines(warn_log, file.path(out_dir, "warnings", paste0(name, ".log")))
  }
  data.frame(
    scenario = name,
    expected = expected,
    status = status,
    pass = isTRUE(pass),
    warn_count = length(warn_log),
    max_warnings = max_warnings,
    nonfinite = isTRUE(nonfinite),
    elapsed_sec = elapsed,
    error = err,
    stringsAsFactors = FALSE
  )
}

stale_plot_args <- c(
  "ci", "deriv", "mean", "plot.view", "intervals", "boot", "bands",
  "plot.errors.method", "plot.errors.type", "plot.errors.alpha",
  "plot.errors.boot.method", "plot.errors.boot.num",
  "plot.errors.boot.nonfixed", "plot.errors.boot.wild",
  "plot.errors.boot.blocklen", "plot.errors.center", "plot.errors.style",
  "plot.errors.bar", "plot.errors.bar.num", "plot.behavior",
  "plot.data.overlay", "plot.rug", "plot.par.mfrow", "plot.bxp",
  "plot.bxp.out", "num.eval", "persp", "xtrim", "xq", "common.scale",
  "display.nomad.progress", "display.warnings"
)

make_stale_case <- function(arg) {
  scenario <- paste0("crs_reject_stale_", gsub("[^[:alnum:]]+", "_", arg))
  substitute(run_case(
    NAME,
    "error",
    "unused plot argument",
    expr = {
      call.args <- list(fit_1d, output = "data")
      call.args[[ARG]] <- TRUE
      do.call(plot, call.args)
    }
  ), list(NAME = scenario, ARG = arg))
}

set.seed(20260611)
n <- 42L
dat <- data.frame(
  x = runif(n),
  x2 = runif(n),
  z = factor(sample(c("a", "b", "c"), n, replace = TRUE))
)
dat$y <- sin(2 * pi * dat$x) + 0.35 * dat$x2 + 0.25 * as.numeric(dat$z) +
  rnorm(n, sd = 0.10)

fit_1d <- crs(
  y ~ x + z,
  data = dat,
  cv = "none",
  kernel = TRUE,
  basis = "additive",
  degree = 2,
  segments = 1,
  lambda = 0.25,
  display.warnings = FALSE,
  display.nomad.progress = FALSE
)

fit_2d <- crs(
  y ~ x + x2,
  data = dat,
  cv = "none",
  kernel = FALSE,
  basis = "additive",
  degree = c(2, 2),
  segments = c(1, 1),
  display.warnings = FALSE,
  display.nomad.progress = FALSE
)

fit_q <- crs(
  y ~ x,
  data = dat,
  cv = "none",
  kernel = FALSE,
  basis = "additive",
  degree = 2,
  segments = 1,
  tau = 0.5,
  display.warnings = FALSE,
  display.nomad.progress = FALSE
)

v <- rnorm(n, sd = 0.08)
w_vec <- runif(n, -1.5, 1.5)
z_vec <- w_vec + v
y_iv <- z_vec^2 + v + rnorm(n, sd = 0.10)
fit_iv <- crsiv(
  y = y_iv,
  z = data.frame(z = z_vec),
  w = data.frame(w = w_vec),
  method = "Tikhonov",
  alpha = 0.1,
  cv = "none",
  basis = "additive",
  deriv = 1,
  display.warnings = FALSE,
  display.nomad.progress = FALSE
)

fit_ivd <- crsivderiv(
  y = y_iv,
  z = data.frame(z = z_vec),
  w = data.frame(w = w_vec),
  iterate.max = 2,
  cv = "none",
  basis = "additive",
  display.warnings = FALSE,
  display.nomad.progress = FALSE
)

fit_clsd <- clsd(
  rnorm(36L),
  degree = 2,
  segments = 1,
  xeval = seq(-2, 2, length.out = 11L),
  display.warnings = FALSE,
  display.nomad.progress = FALSE
)

cases <- list(
  quote(run_case("crs_mean_data_default", "ok",
                 expr = plot(fit_1d, output = "data", neval = neval))),
  quote(run_case("crs_mean_behavior_alias", "ok",
                 expr = plot(fit_1d, behavior = "data", neval = neval))),
  quote(run_case("crs_mean_asymptotic", "ok",
                 expr = plot(fit_1d, output = "data", errors = "asymptotic",
                             band = "pmzsd", neval = neval))),
  quote(run_case("crs_mean_bootstrap_inid", "ok",
                 expr = plot(fit_1d, output = "data", errors = "bootstrap",
                             bootstrap = "inid", B = B, band = "pmzsd",
                             neval = neval,
                             boot_control = np_boot_control(blocklen = 2)))),
  quote(run_case("crs_mean_controls", "ok",
                 expr = plot(fit_1d, output = "data", neval = neval,
                             data_overlay = FALSE, data_rug = TRUE,
                             layout = "auto", common_scale = FALSE,
                             grid_control = np_grid_control(xq = 0.4),
                             render_control = np_render_control(style = "band")))),
  quote(run_case("crs_mean_gradient_data", "ok",
                 expr = plot(fit_1d, output = "data", gradients = TRUE,
                             gradient_order = 1, neval = neval))),
  quote(run_case("crs_mean_gradient_asymptotic", "ok",
                 expr = plot(fit_1d, output = "data", gradients = TRUE,
                             gradient_order = 1, errors = "asymptotic",
                             neval = neval))),
  quote(run_case("crs_mean_render", "ok",
                 expr = run_pdf(plot(fit_1d, neval = neval,
                                     data_overlay = TRUE, main = "CRS fit")))),
  quote(run_case("crs_mean_data_rug_render_clean", "ok",
                 max_warnings = 0L,
                 expr = run_pdf(plot(fit_1d, neval = neval,
                                     data_overlay = TRUE,
                                     data_rug = TRUE)))),
  quote(run_case("crs_surface_data", "ok",
                 expr = plot(fit_2d, output = "data", perspective = TRUE,
                             renderer = "base", neval = 5))),
  quote(run_case("crs_surface_render", "ok",
                 expr = run_pdf(plot(fit_2d, perspective = TRUE,
                                     renderer = "base", neval = 5,
                                     theta = 30, phi = 25,
                                     view = "fixed")))),
  quote(run_case("crs_surface_data_rug_render_clean", "ok",
                 max_warnings = 0L,
                 expr = run_pdf(plot(fit_2d, perspective = TRUE,
                                     renderer = "base", neval = 5,
                                     data_overlay = TRUE,
                                     data_rug = TRUE,
                                     view = "fixed")))),
  quote(run_case("crs_surface_bootstrap_data", "ok",
                 expr = plot(fit_2d, output = "data", perspective = TRUE,
                             errors = "bootstrap", bootstrap = "inid",
                             B = B, neval = 4))),
  quote(run_case("crs_surface_bootstrap_render", "ok",
                 expr = run_pdf(plot(fit_2d, perspective = TRUE,
                                     renderer = "base",
                                     errors = "bootstrap",
                                     bootstrap = "inid", B = B,
                                     neval = 4, data_overlay = TRUE,
                                     data_rug = TRUE,
                                     view = "fixed")))),
  quote(run_case("crs_quantile_data", "ok",
                 expr = plot(fit_q, output = "data", neval = neval))),
  quote(run_case("crs_quantile_gradient", "ok",
                 expr = plot(fit_q, output = "data", gradients = TRUE,
                             gradient_order = 1, neval = neval))),
  quote(run_case("crsiv_phi_data", "ok",
                 expr = plot(fit_iv, output = "data"))),
  quote(run_case("crsiv_derivative_asymptotic", "ok",
                 expr = plot(fit_iv, output = "data", deriv = TRUE,
                             errors = "asymptotic"))),
  quote(run_case("crsiv_render_overlay", "ok",
                 expr = run_pdf(plot(fit_iv, plot.data = TRUE,
                                     main = "CRS IV")))),
  quote(run_case("crsivderiv_default_data", "ok",
                 expr = plot(fit_ivd, output = "data"))),
  quote(run_case("crsivderiv_phi_data", "ok",
                 expr = plot(fit_ivd, output = "data", phi = TRUE))),
  quote(run_case("clsd_density_data", "ok",
                 expr = plot(fit_clsd, output = "data", er = FALSE))),
  quote(run_case("clsd_distribution_data", "ok",
                 expr = plot(fit_clsd, output = "data", distribution = TRUE))),
  quote(run_case("clsd_derivative_render", "ok",
                 expr = run_pdf(plot(fit_clsd, derivative = TRUE,
                                     main = "CLSD derivative")))),
  quote(run_case("crs_reject_ci", "error", "use errors",
                 expr = plot(fit_1d, output = "data", ci = TRUE))),
  quote(run_case("crs_reject_deriv", "error", "use gradients",
                 expr = plot(fit_1d, output = "data", deriv = 1))),
  quote(run_case("crs_reject_mean", "error", "fitted functions",
                 expr = plot(fit_1d, output = "data", mean = TRUE))),
  quote(run_case("crs_reject_plot_view", "error", "NP plot interface",
                 expr = plot(fit_1d, output = "data", plot.view = "fixed"))),
  quote(run_case("crs_reject_boot_without_errors", "error",
                 "bootstrap controls require",
                 expr = plot(fit_1d, output = "data", bootstrap = "inid"))),
  quote(run_case("crs_reject_bad_band_for_asymptotic", "error",
                 "band=\"pmzsd\"",
                 expr = plot(fit_1d, output = "data",
                             errors = "asymptotic", band = "pointwise"))),
  quote(run_case("crs_reject_noninid_bootstrap", "error",
                 "bootstrap=\"inid\"",
                 expr = plot(fit_1d, output = "data",
                             errors = "bootstrap", bootstrap = "wild",
                             B = B))),
  quote(run_case("crs_reject_gradient_bootstrap", "error",
                 "derivative plots",
                 expr = plot(fit_1d, output = "data", gradients = TRUE,
                             errors = "bootstrap", bootstrap = "inid",
                             B = B))),
  quote(run_case("crs_surface_asymptotic_data", "ok",
                 expr = plot(fit_2d, output = "data", perspective = TRUE,
                             errors = "asymptotic", band = "pmzsd",
                             neval = 4))),
  quote(run_case("crs_reject_unknown_arg", "error",
                 "unused plot argument: foo",
                 expr = plot(fit_1d, output = "data",
                             foo = stop("must not evaluate"))))
)

cases <- c(cases, lapply(stale_plot_args, make_stale_case))

rows <- lapply(cases, eval)
res <- do.call(rbind, rows)
res <- res[order(res$scenario), , drop = FALSE]

csv_path <- file.path(out_dir, "crs_plot_stress_results.csv")
write.csv(res, csv_path, row.names = FALSE)

manifest <- c(
  sprintf("repo=%s", repo),
  sprintf("lib=%s", lib),
  sprintf("package_version=%s", as.character(utils::packageVersion("crs"))),
  sprintf("B=%d", B),
  sprintf("neval=%d", neval),
  sprintf("cases=%d", nrow(res)),
  sprintf("failures=%d", sum(!res$pass)),
  sprintf("results_csv=%s", csv_path)
)
writeLines(manifest, file.path(out_dir, "manifest.txt"))

print(res[, c("scenario", "expected", "status", "pass", "warn_count",
              "nonfinite", "elapsed_sec", "error")], row.names = FALSE)
cat(sprintf("results_csv=%s\n", csv_path))
cat(sprintf("manifest=%s\n", file.path(out_dir, "manifest.txt")))

if (!all(res$pass)) quit(save = "no", status = 1L)
quit(save = "no", status = 0L)
