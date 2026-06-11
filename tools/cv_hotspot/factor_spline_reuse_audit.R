args <- commandArgs(trailingOnly = TRUE)
label <- if (length(args) >= 1L) args[[1L]] else "candidate"
out_dir <- if (length(args) >= 2L) args[[2L]] else file.path(tempdir(), "factor_spline_reuse_audit")
n <- if (length(args) >= 3L) as.integer(args[[3L]]) else 4000L
reps <- if (length(args) >= 4L) as.integer(args[[4L]]) else 120L

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(crs))

set.seed(20260611L)
x <- cbind(runif(n), runif(n))
z <- data.frame(
  f1 = factor(sample(letters[1:4], n, replace = TRUE)),
  f2 = factor(sample(LETTERS[1:5], n, replace = TRUE))
)
f1_eff <- c(a = -0.2, b = 0.0, c = 0.25, d = 0.45)[as.character(z$f1)]
f2_eff <- c(A = -0.15, B = 0.1, C = 0.0, D = 0.2, E = -0.05)[as.character(z$f2)]
signal <- sin(2 * pi * x[, 1L]) + cos(2 * pi * x[, 2L]) + f1_eff + f2_eff
y <- signal + rnorm(n, sd = 0.4 * sd(signal))
w <- 0.35 + runif(n)

K <- matrix(c(4, 2, 3, 2), ncol = 2L, byrow = TRUE)
I <- c(1L, 1L)
rows <- seq_len(n)

time_loop <- function(expr_fun, reps) {
  gc(FALSE)
  invisible(expr_fun())
  elapsed <- system.time({
    for (i in seq_len(reps)) invisible(expr_fun())
  })[["elapsed"]]
  elapsed / reps
}

cv_from_fit <- function(fit, cv.func, y, weights) {
  eps <- fit$residuals.rows
  htt <- fit$hat.rows
  nobs <- length(y)
  if (!is.null(weights)) eps <- eps * sqrt(weights)
  if (is.null(weights)) {
    if (cv.func == "cv.ls") {
      mean((eps / (1 - htt))^2)
    } else if (cv.func == "cv.gcv") {
      mean(eps^2 / (1 - mean(htt))^2)
    } else {
      traceH <- sum(htt)
      penalty <- ((1 + traceH / nobs) / (1 - (traceH + 2) / nobs))
      if (penalty < 0) crs:::resolve_cv_maxPenalty(NULL, y, weights = weights, cv.func = cv.func) else log(mean(eps^2)) + penalty
    }
  } else {
    if (cv.func == "cv.ls") {
      mean(eps^2 / (1 - htt)^2)
    } else if (cv.func == "cv.gcv") {
      mean(eps^2 / (1 - mean(htt))^2)
    } else {
      traceH <- sum(htt)
      penalty <- ((1 + traceH / nobs) / (1 - (traceH + 2) / nobs))
      if (penalty < 0) crs:::resolve_cv_maxPenalty(NULL, y, weights = weights, cv.func = cv.func) else log(mean(eps^2)) + penalty
    }
  }
}

run_one <- function(basis, weighted, cv.func) {
  weights <- if (weighted) w else NULL
  P <- crs:::prod.spline(
    x = x,
    z = z,
    K = K,
    I = I,
    knots = "quantiles",
    basis = basis,
    display.warnings = FALSE
  )
  X <- if (basis %in% c("additive", "glp")) cbind(1, P) else P

  full_cv <- crs:::cv.factor.spline(
    x = x,
    y = y,
    z = z,
    K = K,
    I = I,
    knots = "quantiles",
    basis = basis,
    cv.func = cv.func,
    weights = weights,
    singular.ok = FALSE,
    display.warnings = FALSE,
    use.ridge = FALSE,
    use.gram.cv = TRUE
  )

  fit <- crs:::.crs_weighted_ls_cv_rows(
    X = X,
    y = y,
    weights = weights,
    rows = rows,
    ridge.lambda = NULL,
    rcond.min = 1e-8,
    allow.fallback = TRUE,
    use.svd.fallback = TRUE
  )
  reconstructed_cv <- cv_from_fit(fit, cv.func, y, weights)

  full_time <- time_loop(function() {
    crs:::cv.factor.spline(
      x = x,
      y = y,
      z = z,
      K = K,
      I = I,
      knots = "quantiles",
      basis = basis,
      cv.func = cv.func,
      weights = weights,
      singular.ok = FALSE,
      display.warnings = FALSE,
      use.ridge = FALSE,
      use.gram.cv = TRUE
    )
  }, reps)

  basis_time <- time_loop(function() {
    crs:::prod.spline(
      x = x,
      z = z,
      K = K,
      I = I,
      knots = "quantiles",
      basis = basis,
      display.warnings = FALSE
    )
  }, reps)

  ls_time <- time_loop(function() {
    crs:::.crs_weighted_ls_cv_rows(
      X = X,
      y = y,
      weights = weights,
      rows = rows,
      ridge.lambda = NULL,
      rcond.min = 1e-8,
      allow.fallback = TRUE,
      use.svd.fallback = TRUE
    )
  }, reps)

  factor_mm_time <- time_loop(function() {
    for (j in seq_along(z)) invisible(model.matrix(~z[, j])[, -1L, drop = FALSE])
  }, reps)

  data.frame(
    label = label,
    n = n,
    reps = reps,
    basis = basis,
    weighted = weighted,
    cv.func = cv.func,
    ncol.P = ncol(P),
    ncol.X = ncol(X),
    full_cv = as.numeric(full_cv),
    reconstructed_cv = as.numeric(reconstructed_cv),
    abs_delta = abs(as.numeric(full_cv) - as.numeric(reconstructed_cv)),
    full_time = full_time,
    basis_time = basis_time,
    ls_time = ls_time,
    factor_mm_time = factor_mm_time,
    basis_share = basis_time / full_time,
    ls_share = ls_time / full_time,
    factor_mm_share = factor_mm_time / full_time,
    stringsAsFactors = FALSE
  )
}

grid <- expand.grid(
  basis = c("additive", "tensor", "glp"),
  weighted = c(FALSE, TRUE),
  cv.func = c("cv.ls", "cv.gcv", "cv.aic"),
  stringsAsFactors = FALSE
)

results <- do.call(rbind, Map(
  run_one,
  basis = grid$basis,
  weighted = grid$weighted,
  cv.func = grid$cv.func
))

csv_path <- file.path(out_dir, paste0("factor_spline_reuse_audit_", label, ".csv"))
write.csv(results, csv_path, row.names = FALSE)

summary_rows <- aggregate(
  cbind(full_time, basis_time, ls_time, basis_share, ls_share, factor_mm_share) ~ basis + weighted,
  data = results,
  FUN = median
)
summary_path <- file.path(out_dir, paste0("factor_spline_reuse_audit_summary_", label, ".csv"))
write.csv(summary_rows, summary_path, row.names = FALSE)

md_path <- file.path(out_dir, paste0("factor_spline_reuse_audit_", label, ".md"))
con <- file(md_path, open = "wt")
on.exit(close(con), add = TRUE)
writeLines(c(
  "# Factor-Spline CV Reuse Audit",
  "",
  paste0("Label: `", label, "`"),
  paste0("n: `", n, "`"),
  paste0("reps per timing cell: `", reps, "`"),
  "",
  "## Numerical Check",
  "",
  paste0("Maximum absolute full-vs-reconstructed CV delta: `", signif(max(results$abs_delta), 8), "`"),
  "",
  "## Median Timing Shares",
  "",
  paste(capture.output(print(summary_rows, row.names = FALSE)), collapse = "\n"),
  "",
  "## Files",
  "",
  paste0("- Raw CSV: `", csv_path, "`"),
  paste0("- Summary CSV: `", summary_path, "`")
), con)

cat("wrote", csv_path, "\n")
cat("wrote", summary_path, "\n")
cat("wrote", md_path, "\n")
