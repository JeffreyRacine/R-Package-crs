#!/usr/bin/env Rscript

suppressWarnings(suppressPackageStartupMessages(library(crs)))

args <- commandArgs(trailingOnly = TRUE)
out_prefix <- if (length(args) >= 1) args[[1]] else "/tmp/crs_nomad_solver_space"
n <- if (length(args) >= 2) as.integer(args[[2]]) else 180L

mk_data <- function(scn, n, seed) {
  set.seed(seed)
  if (scn == "simple") {
    x <- runif(n)
    z <- factor(sample(c("A", "B"), n, replace = TRUE))
    y <- sin(2 * pi * x) + as.numeric(z == "B") + rnorm(n, sd = 0.2)
    return(list(
      frkr_xz = data.frame(x = x, z = z),
      frkr_y = y,
      np_df = data.frame(y = y, x = x, z = z),
      np_formula = y ~ x + z,
      seg = c(2),
      dmax = 4,
      dmin = 0
    ))
  }

  x1 <- runif(n)
  x2 <- runif(n)
  z <- factor(sample(c("A", "B", "C"), n, replace = TRUE, prob = c(0.35, 0.35, 0.30)))
  y <- sin(2 * pi * x1) + cos(2 * pi * x2) + 0.5 * x1 * x2 +
    as.numeric(z == "B") * 0.4 - as.numeric(z == "C") * 0.25 + rnorm(n, sd = 0.2)
  list(
    frkr_xz = data.frame(x1 = x1, x2 = x2, z = z),
    frkr_y = y,
    np_df = data.frame(y = y, x1 = x1, x2 = x2, z = z),
    np_formula = y ~ x1 + x2 + z,
    seg = c(2, 2),
    dmax = 4,
    dmin = 0
  )
}

run_case <- function(case_id, d, strategy) {
  nmulti <- as.integer(strategy$nmulti)
  maxbb <- as.integer(strategy$max_bb_eval)
  opts <- strategy$opts
  tm <- system.time({
    out <- switch(
      case_id,
      frscvNOMAD = frscvNOMAD(
        xz = d$frkr_xz,
        y = d$frkr_y,
        complexity = "degree",
        segments = d$seg,
        degree.max = d$dmax,
        degree.min = d$dmin,
        nmulti = nmulti,
        max.bb.eval = maxbb,
        display.warnings = FALSE,
        display.nomad.progress = FALSE,
        opts = opts
      ),
      krscvNOMAD = krscvNOMAD(
        xz = d$frkr_xz,
        y = d$frkr_y,
        complexity = "degree",
        segments = d$seg,
        degree.max = d$dmax,
        degree.min = d$dmin,
        nmulti = nmulti,
        max.bb.eval = maxbb,
        display.warnings = FALSE,
        display.nomad.progress = FALSE,
        opts = opts
      ),
      npglpreg = npglpreg(
        formula = d$np_formula,
        data = d$np_df,
        cv = "degree-bandwidth",
        nmulti = nmulti,
        max.bb.eval = maxbb,
        degree.max = d$dmax,
        degree.min = d$dmin,
        display.warnings = FALSE,
        display.nomad.progress = FALSE,
        opts = opts
      )
    )
  })[["elapsed"]]

  obj <- switch(
    case_id,
    frscvNOMAD = out$cv.objc[1],
    krscvNOMAD = out$cv.objc[1],
    npglpreg = out$fv[1]
  )

  c(elapsed = as.numeric(tm), objective = as.numeric(obj))
}

strategies <- list(
  frscvNOMAD = list(
    fr_current = list(nmulti = 2L, max_bb_eval = 100L, opts = list()),
    fr_noquad_lexi = list(nmulti = 2L, max_bb_eval = 100L, opts = list("QUAD_MODEL_SEARCH" = "no", "EVAL_QUEUE_SORT" = "LEXICOGRAPHICAL")),
    fr_noquad_random = list(nmulti = 2L, max_bb_eval = 100L, opts = list("QUAD_MODEL_SEARCH" = "no", "EVAL_QUEUE_SORT" = "RANDOM")),
    fr_cs_fast = list(nmulti = 2L, max_bb_eval = 100L, opts = list("CS_OPTIMIZATION" = "yes")),
    fr_cs_eval200 = list(nmulti = 2L, max_bb_eval = 200L, opts = list("CS_OPTIMIZATION" = "yes")),
    fr_cs_nmulti4 = list(nmulti = 4L, max_bb_eval = 100L, opts = list("CS_OPTIMIZATION" = "yes")),
    fr_cs_eval200_nmulti4 = list(nmulti = 4L, max_bb_eval = 200L, opts = list("CS_OPTIMIZATION" = "yes")),
    fr_random_opt = list(nmulti = 2L, max_bb_eval = 100L, opts = list("RANDOM_ALGO_OPTIMIZATION" = "yes")),
    fr_nm_opt = list(nmulti = 2L, max_bb_eval = 100L, opts = list("NM_OPTIMIZATION" = "yes")),
    fr_qp_opt = list(nmulti = 2L, max_bb_eval = 100L, opts = list("QP_OPTIMIZATION" = "yes"))
  ),
  krscvNOMAD = list(
    kr_current = list(nmulti = 2L, max_bb_eval = 100L, opts = list()),
    kr_noquad_lexi = list(nmulti = 2L, max_bb_eval = 100L, opts = list("QUAD_MODEL_SEARCH" = "no", "EVAL_QUEUE_SORT" = "LEXICOGRAPHICAL")),
    kr_noquad_random = list(nmulti = 2L, max_bb_eval = 100L, opts = list("QUAD_MODEL_SEARCH" = "no", "EVAL_QUEUE_SORT" = "RANDOM")),
    kr_cs_fast = list(nmulti = 2L, max_bb_eval = 100L, opts = list("CS_OPTIMIZATION" = "yes")),
    kr_cs_eval200 = list(nmulti = 2L, max_bb_eval = 200L, opts = list("CS_OPTIMIZATION" = "yes")),
    kr_cs_nmulti4 = list(nmulti = 4L, max_bb_eval = 100L, opts = list("CS_OPTIMIZATION" = "yes")),
    kr_cs_eval200_nmulti4 = list(nmulti = 4L, max_bb_eval = 200L, opts = list("CS_OPTIMIZATION" = "yes")),
    kr_random_opt = list(nmulti = 2L, max_bb_eval = 100L, opts = list("RANDOM_ALGO_OPTIMIZATION" = "yes")),
    kr_nm_opt = list(nmulti = 2L, max_bb_eval = 100L, opts = list("NM_OPTIMIZATION" = "yes")),
    kr_qp_opt = list(nmulti = 2L, max_bb_eval = 100L, opts = list("QP_OPTIMIZATION" = "yes"))
  ),
  npglpreg = list(
    np_current = list(nmulti = 2L, max_bb_eval = 100L, opts = list()),
    np_quadbox1 = list(nmulti = 2L, max_bb_eval = 100L, opts = list("QUAD_MODEL_SEARCH_BOX_FACTOR" = "1.0", "QUAD_MODEL_BOX_FACTOR" = "1.0")),
    np_noquad = list(nmulti = 2L, max_bb_eval = 100L, opts = list("QUAD_MODEL_SEARCH" = "no", "EVAL_QUEUE_SORT" = "DIR_LAST_SUCCESS")),
    np_noquad_eval200 = list(nmulti = 2L, max_bb_eval = 200L, opts = list("QUAD_MODEL_SEARCH" = "no", "EVAL_QUEUE_SORT" = "DIR_LAST_SUCCESS")),
    np_noquad_nmulti4 = list(nmulti = 4L, max_bb_eval = 100L, opts = list("QUAD_MODEL_SEARCH" = "no", "EVAL_QUEUE_SORT" = "DIR_LAST_SUCCESS")),
    np_noquad_eval200_nmulti4 = list(nmulti = 4L, max_bb_eval = 200L, opts = list("QUAD_MODEL_SEARCH" = "no", "EVAL_QUEUE_SORT" = "DIR_LAST_SUCCESS")),
    np_cs_fast = list(nmulti = 2L, max_bb_eval = 100L, opts = list("CS_OPTIMIZATION" = "yes")),
    np_cs_eval200 = list(nmulti = 2L, max_bb_eval = 200L, opts = list("CS_OPTIMIZATION" = "yes")),
    np_nm_opt = list(nmulti = 2L, max_bb_eval = 100L, opts = list("NM_OPTIMIZATION" = "yes")),
    np_qp_opt = list(nmulti = 2L, max_bb_eval = 100L, opts = list("QP_OPTIMIZATION" = "yes")),
    np_random_opt = list(nmulti = 2L, max_bb_eval = 100L, opts = list("RANDOM_ALGO_OPTIMIZATION" = "yes"))
  )
)

baselines <- c(frscvNOMAD = "fr_current", krscvNOMAD = "kr_current", npglpreg = "np_current")
tols <- c(frscvNOMAD = 1e-10, krscvNOMAD = 1e-6, npglpreg = 5e-4)

rows <- list()
idx <- 1L
for (scn in c("simple", "hard")) {
  cat("SCENARIO", scn, "\n")
  for (spol in c("fixed_seed", "varying_seed")) {
    seeds <- if (spol == "fixed_seed") rep(42L, 5L) else c(1001L, 1002L, 1003L, 1004L, 1005L)
    for (sd in seeds) {
      d <- mk_data(scn, n = n, seed = sd)
      for (case_id in names(strategies)) {
        for (sname in names(strategies[[case_id]])) {
          st <- strategies[[case_id]][[sname]]
          status <- "ok"
          err <- ""
          el <- NA_real_
          obj <- NA_real_
          rr <- tryCatch(
            run_case(case_id, d, st),
            error = function(e) {
              status <<- "error"
              err <<- conditionMessage(e)
              NULL
            }
          )
          if (!is.null(rr)) {
            el <- rr[["elapsed"]]
            obj <- rr[["objective"]]
          }
          rows[[idx]] <- data.frame(
            scenario = scn,
            seed_policy = spol,
            seed = as.integer(sd),
            case = case_id,
            profile = sname,
            nmulti = as.integer(st$nmulti),
            max_bb_eval = as.integer(st$max_bb_eval),
            status = status,
            elapsed_sec = el,
            objective = obj,
            error = err,
            stringsAsFactors = FALSE
          )
          idx <- idx + 1L
        }
      }
    }
  }
}

raw <- do.call(rbind, rows)
raw$replicate <- ave(
  seq_len(nrow(raw)),
  interaction(raw$scenario, raw$case, raw$profile, raw$seed_policy, raw$seed, drop = TRUE),
  FUN = seq_along
)

ok <- raw[raw$status == "ok" & is.finite(raw$objective) & is.finite(raw$elapsed_sec), ]
cmp_rows <- list()
cidx <- 1L
for (scn in unique(ok$scenario)) {
  for (case_id in unique(ok$case)) {
    base_name <- baselines[[case_id]]
    tol <- tols[[case_id]]
    for (spol in unique(ok$seed_policy)) {
      base <- ok[
        ok$scenario == scn & ok$case == case_id & ok$profile == base_name & ok$seed_policy == spol,
        c("seed", "replicate", "elapsed_sec", "objective")
      ]
      names(base)[3:4] <- c("elapsed_base", "objective_base")
      if (nrow(base) == 0) next

      for (pname in unique(ok$profile[ok$scenario == scn & ok$case == case_id])) {
        cand <- ok[
          ok$scenario == scn & ok$case == case_id & ok$profile == pname & ok$seed_policy == spol,
          c("seed", "replicate", "elapsed_sec", "objective")
        ]
        names(cand)[3:4] <- c("elapsed_cand", "objective_cand")
        if (nrow(cand) == 0) next

        m <- merge(base, cand, by = c("seed", "replicate"), all = FALSE)
        if (nrow(m) == 0) next

        dp <- 100 * (m$elapsed_cand - m$elapsed_base) / pmax(m$elapsed_base, .Machine$double.eps)
        od <- m$objective_cand - m$objective_base

        cmp_rows[[cidx]] <- data.frame(
          scenario = scn,
          case = case_id,
          seed_policy = spol,
          profile = pname,
          paired_runs = nrow(m),
          mean_elapsed_pct_change = mean(dp),
          median_elapsed_pct_change = median(dp),
          mean_objective_diff = mean(od),
          median_objective_diff = median(od),
          worsen_count = sum(od > tol, na.rm = TRUE),
          improve_count = sum(od < -tol, na.rm = TRUE),
          tie_count = sum(abs(od) <= tol, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
        cidx <- cidx + 1L
      }
    }
  }
}
cmp <- do.call(rbind, cmp_rows)

agg <- aggregate(cbind(worsen_count, improve_count) ~ case + profile, data = cmp, FUN = sum)
coverage <- aggregate(mean_elapsed_pct_change ~ case + profile, data = cmp, FUN = length)
names(coverage)[3] <- "cells_present"
mme <- aggregate(mean_elapsed_pct_change ~ case + profile, data = cmp, FUN = mean)
med <- aggregate(mean_elapsed_pct_change ~ case + profile, data = cmp, FUN = median)
worst <- aggregate(mean_elapsed_pct_change ~ case + profile, data = cmp, FUN = max)
best <- aggregate(mean_elapsed_pct_change ~ case + profile, data = cmp, FUN = min)
agg <- merge(agg, coverage, by = c("case", "profile"), all.x = TRUE)
agg <- merge(agg, mme, by = c("case", "profile"), all.x = TRUE)
names(agg)[names(agg) == "mean_elapsed_pct_change"] <- "mean_elapsed_pct_change_overall"
agg$median_elapsed_pct_change_overall <- med$mean_elapsed_pct_change
agg$worst_case_elapsed_pct_change <- worst$mean_elapsed_pct_change
agg$best_case_elapsed_pct_change <- best$mean_elapsed_pct_change

err <- raw[raw$status != "ok", ]
if (nrow(err) > 0) {
  er <- aggregate(status ~ case + profile, data = err, FUN = length)
  names(er)[3] <- "error_runs"
  agg <- merge(agg, er, by = c("case", "profile"), all.x = TRUE)
} else {
  agg$error_runs <- 0L
}
agg$error_runs[is.na(agg$error_runs)] <- 0L

agg <- agg[order(agg$case, agg$worsen_count, agg$error_runs, agg$mean_elapsed_pct_change_overall), ]
strict_safe <- agg[agg$cells_present == 4 & agg$worsen_count == 0 & agg$error_runs == 0, ]
strict_safe <- strict_safe[order(strict_safe$case, strict_safe$mean_elapsed_pct_change_overall), ]

raw_path <- paste0(out_prefix, "_raw.csv")
cmp_path <- paste0(out_prefix, "_compare.csv")
agg_path <- paste0(out_prefix, "_agg.csv")
safe_path <- paste0(out_prefix, "_strict_safe.csv")

write.csv(raw, raw_path, row.names = FALSE)
write.csv(cmp, cmp_path, row.names = FALSE)
write.csv(agg, agg_path, row.names = FALSE)
write.csv(strict_safe, safe_path, row.names = FALSE)

cat("WROTE", raw_path, "\n")
cat("WROTE", cmp_path, "\n")
cat("WROTE", agg_path, "\n")
cat("WROTE", safe_path, "\n")
cat("ERROR_RUNS", nrow(err), "\n")
