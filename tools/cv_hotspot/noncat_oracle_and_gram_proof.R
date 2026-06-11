#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if(length(args)) args[[1L]] else
  "/Users/jracine/Development/tmp/crs_cv_hotspot_engineering_20260611/oracles/noncat_gram"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(crs))

timer <- function(expr) {
  gc()
  t0 <- proc.time()
  value <- force(expr)
  elapsed <- unname((proc.time() - t0)[["elapsed"]])
  list(value = value, elapsed = elapsed)
}

make_data <- function(n = 160L, p = 2L, seed = 1L, rank_stress = FALSE) {
  set.seed(seed)
  x <- matrix(runif(n * p), ncol = p)
  colnames(x) <- paste0("x", seq_len(p))
  if(rank_stress && p >= 2L) {
    x[, p] <- x[, 1L] + rnorm(n, sd = 1e-9)
  }
  signal <- sin(2 * pi * x[, 1L])
  if(p >= 2L) signal <- signal + 0.5 * cos(2 * pi * x[, 2L])
  y <- signal + rnorm(n, sd = 0.2)
  weights <- seq(0.6, 1.4, length.out = n)
  list(x = x, y = y, weights = weights)
}

cv_current <- function(dat, basis, cv_func, weighted, K, use_gram = TRUE) {
  crs:::cv.kernel.spline(
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
    weights = if(weighted) dat$weights else NULL,
    tau = NULL,
    singular.ok = FALSE,
    display.warnings = FALSE,
    use.ridge = FALSE,
    smooth.penalty = TRUE,
    penalty.scale = 1000,
    use.gram.cv = use_gram,
    gram.rcond.min = 1e-10,
    record.gram.stats = TRUE
  )
}

cv_from_payload <- function(epsilon, htt, cv_func, y, weights = NULL) {
  if(!is.null(weights)) epsilon <- epsilon * sqrt(weights)
  n <- length(y)
  if(cv_func == "cv.ls") {
    mean(epsilon^2 / (1 - htt)^2)
  } else if(cv_func == "cv.gcv") {
    mean(epsilon^2 / (1 - mean(htt))^2)
  } else if(cv_func == "cv.aic") {
    traceH <- sum(htt)
    sigmasq <- mean(epsilon^2)
    penalty <- ((1 + traceH / n) / (1 - (traceH + 2) / n))
    if(penalty < 0) Inf else log(sigmasq) + penalty
  } else {
    stop("unknown cv_func")
  }
}

noncat_payloads <- function(dat, basis, K, weighted) {
  P <- crs:::prod.spline(x = dat$x, K = K, knots = "quantiles",
                         basis = basis, display.warnings = FALSE)
  X <- if(basis %in% c("additive", "glp")) cbind(1, P) else P
  weights <- if(weighted) dat$weights else NULL

  fit <- crs:::.crs_weighted_ls_cv_rows(
    X = X,
    y = dat$y,
    weights = weights,
    rows = seq_len(NROW(X)),
    rcond.min = 1e-10,
    ridge.lambda = NULL,
    allow.fallback = TRUE,
    use.svd.fallback = TRUE
  )

  core <- crs:::.crs_weighted_ls_core(
    X = X,
    y = dat$y,
    weights = weights,
    rcond.min = 1e-10,
    ridge.lambda = NULL,
    allow.fallback = TRUE,
    use.svd.fallback = TRUE
  )
  epsilon_core <- dat$y - drop(X %*% core$coefficients)

  list(
    P = P,
    X = X,
    weighted = weighted,
    weights = weights,
    fit = fit,
    epsilon_core = epsilon_core,
    hat_legacy = pmin(hat(P), 1 - .Machine$double.eps),
    hat_design = pmin(fit$hat.rows, 1 - .Machine$double.eps)
  )
}

cases <- expand.grid(
  basis = c("additive", "tensor", "glp"),
  cv_func = c("cv.ls", "cv.gcv", "cv.aic"),
  weighted = c(FALSE, TRUE),
  rank_stress = c(FALSE, TRUE),
  stringsAsFactors = FALSE
)

K_regular <- matrix(c(3, 2, 2, 2), ncol = 2L, byrow = TRUE)
K_rank <- matrix(c(4, 2, 4, 2), ncol = 2L, byrow = TRUE)

results <- vector("list", nrow(cases))
for(i in seq_len(nrow(cases))) {
  row <- cases[i, ]
  dat <- make_data(n = if(row$rank_stress) 120L else 180L,
                   p = 2L,
                   seed = 81000L + i,
                   rank_stress = row$rank_stress)
  K <- if(row$rank_stress) K_rank else K_regular
  current <- timer(cv_current(dat, row$basis, row$cv_func, row$weighted, K))
  payload <- noncat_payloads(dat, row$basis, K, row$weighted)
  legacy_cv <- cv_from_payload(payload$epsilon_core, payload$hat_legacy,
                               row$cv_func, dat$y, payload$weights)
  design_cv <- cv_from_payload(payload$epsilon_core, payload$hat_design,
                               row$cv_func, dat$y, payload$weights)
  results[[i]] <- cbind(
    cases[i, ],
    data.frame(
      n = NROW(dat$x),
      p_design = NCOL(payload$X),
      current = as.numeric(current$value),
      legacy_reconstructed = legacy_cv,
      design_reconstructed = design_cv,
      abs_current_legacy = abs(as.numeric(current$value) - legacy_cv),
      abs_current_design = abs(as.numeric(current$value) - design_cv),
      max_abs_hat_delta = max(abs(payload$hat_legacy - payload$hat_design)),
      trace_legacy = sum(payload$hat_legacy),
      trace_design = sum(payload$hat_design),
      method = payload$fit$method,
      status = payload$fit$status,
      rcond = payload$fit$rcond,
      elapsed_current = current$elapsed
    )
  )
}
summary_df <- do.call(rbind, results)

route_map <- data.frame(
  route = c("cv.kernel.spline z=NULL mean",
            "cv.kernel.spline z!=NULL mean",
            "cv.factor.spline mean",
            "cv.kernel.spline quantile tau",
            "cv.factor.spline quantile tau"),
  current_solver = c("LS core for coefficients; legacy hat(P) leverage",
                     "weighted-LS rows/system; cell cache when eligible",
                     "weighted-LS rows",
                     "rq plus lm hat for leverage",
                     "rq plus lm hat for leverage"),
  tranche_status = c("Tranche B/C target",
                     "Already Gram/cell-cache candidate; later kernel/native decisions",
                     "Tranche D structural audit",
                     "Out of scope for first mean-regression performance tranche",
                     "Out of scope for first mean-regression performance tranche"),
  stringsAsFactors = FALSE
)

saveRDS(list(route_map = route_map,
             cases = cases,
             summary = summary_df,
             session = sessionInfo()),
        file.path(out_dir, "noncat_oracle_and_gram_proof.rds"))
write.csv(summary_df, file.path(out_dir, "noncat_oracle_and_gram_proof.csv"),
          row.names = FALSE)
writeLines(c(
  "# Non-Categorical Gram Oracle And Proof",
  "",
  paste("Date:", as.character(Sys.time())),
  "",
  "## Route Map",
  "",
  paste(capture.output(print(route_map, row.names = FALSE)), collapse = "\n"),
  "",
  "## Summary",
  "",
  paste(capture.output(print(summary_df, row.names = FALSE, digits = 12)),
        collapse = "\n"),
  "",
  "## Max Deltas",
  "",
  sprintf("max abs current-vs-legacy: %.17g",
          max(summary_df$abs_current_legacy, na.rm = TRUE)),
  sprintf("max abs current-vs-design: %.17g",
          max(summary_df$abs_current_design, na.rm = TRUE)),
  sprintf("max abs legacy/design hat delta: %.17g",
          max(summary_df$max_abs_hat_delta, na.rm = TRUE))
), file.path(out_dir, "NONCAT_GRAM_ORACLE_AND_PROOF.md"))

print(summary_df, row.names = FALSE, digits = 12)
cat("\nWrote artifacts to", out_dir, "\n")
