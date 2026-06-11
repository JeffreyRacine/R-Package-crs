#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
label <- if(length(args) >= 1L) args[[1L]] else "run"
out_dir <- if(length(args) >= 2L) args[[2L]] else
  "/Users/jracine/Development/tmp/crs_cv_hotspot_engineering_20260611/timing/noncat_fixed_eval"
reps <- if(length(args) >= 3L) as.integer(args[[3L]]) else 160L
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(crs))

make_data <- function(n = 4000L, seed = 91001L) {
  set.seed(seed)
  x <- cbind(x1 = runif(n), x2 = runif(n))
  signal <- sin(2 * pi * x[, 1L]) + 0.5 * cos(2 * pi * x[, 2L])
  y <- signal + rnorm(n, sd = 0.25)
  weights <- seq(0.6, 1.4, length.out = n)
  list(x = x, y = y, weights = weights)
}

run_one <- function(dat, basis, cv_func, weighted, reps) {
  K <- matrix(c(4, 2, 3, 2), ncol = 2L, byrow = TRUE)
  weights <- if(weighted) dat$weights else NULL
  values <- numeric(reps)
  gc()
  t0 <- proc.time()
  for(i in seq_len(reps)) {
    values[[i]] <- crs:::cv.kernel.spline(
      x = dat$x,
      y = dat$y,
      z = NULL,
      K = K,
      lambda = NULL,
      z.unique = NULL,
      ind = NULL,
      ind.vals = NULL,
      ind.list = NULL,
      nrow.z.unique = 0L,
      is.ordered.z = FALSE,
      knots = "quantiles",
      basis = basis,
      cv.func = cv_func,
      weights = weights,
      tau = NULL,
      singular.ok = FALSE,
      display.warnings = FALSE,
      use.ridge = FALSE,
      smooth.penalty = TRUE,
      penalty.scale = 1000,
      use.gram.cv = TRUE,
      gram.rcond.min = 1e-10,
      record.gram.stats = FALSE
    )
  }
  elapsed <- unname((proc.time() - t0)[["elapsed"]])
  data.frame(
    label = label,
    basis = basis,
    cv_func = cv_func,
    weighted = weighted,
    reps = reps,
    value_first = values[[1L]],
    value_last = values[[reps]],
    max_abs_value_delta_internal = max(abs(values - values[[1L]])),
    elapsed = elapsed,
    per_eval = elapsed / reps,
    stringsAsFactors = FALSE
  )
}

dat <- make_data()
cases <- expand.grid(
  basis = c("additive", "tensor", "glp"),
  cv_func = c("cv.ls", "cv.gcv", "cv.aic"),
  weighted = c(FALSE, TRUE),
  stringsAsFactors = FALSE
)
res <- do.call(rbind, lapply(seq_len(nrow(cases)), function(i) {
  row <- cases[i, ]
  run_one(dat, row$basis, row$cv_func, row$weighted, reps)
}))
out <- file.path(out_dir, paste0(label, "_noncat_fixed_eval_timing.csv"))
write.csv(res, out, row.names = FALSE)
print(res, row.names = FALSE, digits = 12)
cat("\nWrote", out, "\n")
