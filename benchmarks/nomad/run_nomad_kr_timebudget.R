#!/usr/bin/env Rscript

suppressWarnings(suppressPackageStartupMessages(library(crs)))

args <- commandArgs(trailingOnly = TRUE)
out_prefix <- if (length(args) >= 1) args[[1]] else "/tmp/crs_nomad_kr_timebudget"
n <- if (length(args) >= 2) as.integer(args[[2]]) else 300L

mk_data <- function(scn, n, seed) {
  set.seed(seed)
  if (scn == "simple") {
    x <- runif(n)
    z <- factor(sample(c("A", "B"), n, replace = TRUE))
    y <- sin(2 * pi * x) + as.numeric(z == "B") + rnorm(n, sd = 0.2)
    return(list(
      xz = data.frame(x = x, z = z),
      y = y,
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
    xz = data.frame(x1 = x1, x2 = x2, z = z),
    y = y,
    seg = c(2, 2),
    dmax = 4,
    dmin = 0
  )
}

run_kr <- function(d, cfg) {
  tm <- system.time({
    out <- krscvNOMAD(
      xz = d$xz,
      y = d$y,
      complexity = "degree",
      segments = d$seg,
      degree.max = d$dmax,
      degree.min = d$dmin,
      nmulti = as.integer(cfg$nmulti),
      max.bb.eval = as.integer(cfg$max_bb_eval),
      opts = cfg$opts,
      display.warnings = FALSE,
      display.nomad.progress = FALSE
    )
  })[["elapsed"]]

  list(
    elapsed = as.numeric(tm),
    objective = as.numeric(out$cv.objc[1]),
    degree = paste(out$degree, collapse = ";"),
    lambda = paste(signif(out$lambda, 12), collapse = ";")
  )
}

profiles <- list(
  kr_current = list(nmulti = 2L, max_bb_eval = 100L, opts = list()),
  kr_nmtrial40 = list(nmulti = 2L, max_bb_eval = 100L, opts = list("NM_SEARCH_MAX_TRIAL_PTS_NFACTOR" = 40L)),
  kr_dirnp1neg = list(nmulti = 2L, max_bb_eval = 100L, opts = list("DIRECTION_TYPE" = "ORTHO N+1 NEG")),
  kr_dir2n = list(nmulti = 2L, max_bb_eval = 100L, opts = list("DIRECTION_TYPE" = "ORTHO 2N")),
  kr_simple_specoff = list(nmulti = 2L, max_bb_eval = 100L, opts = list("SIMPLE_LINE_SEARCH" = "yes", "SPECULATIVE_SEARCH" = "no")),
  kr_eval80 = list(nmulti = 2L, max_bb_eval = 80L, opts = list()),
  kr_dirnp1neg_nmulti4 = list(nmulti = 4L, max_bb_eval = 100L, opts = list("DIRECTION_TYPE" = "ORTHO N+1 NEG")),
  kr_dir2n_nmulti4 = list(nmulti = 4L, max_bb_eval = 100L, opts = list("DIRECTION_TYPE" = "ORTHO 2N")),
  kr_dirnp1neg_eval140 = list(nmulti = 2L, max_bb_eval = 140L, opts = list("DIRECTION_TYPE" = "ORTHO N+1 NEG")),
  kr_dir2n_eval140 = list(nmulti = 2L, max_bb_eval = 140L, opts = list("DIRECTION_TYPE" = "ORTHO 2N"))
)

rows <- list()
idx <- 1L
for (scn in c("simple", "hard")) {
  cat("SCENARIO", scn, "\n")
  for (spol in c("fixed_seed", "varying_seed")) {
    seeds <- if (spol == "fixed_seed") rep(42L, 5L) else c(1001L, 1002L, 1003L, 1004L, 1005L)
    for (sd in seeds) {
      d <- mk_data(scn, n = n, seed = sd)
      for (pname in names(profiles)) {
        cfg <- profiles[[pname]]
        status <- "ok"
        err <- ""
        el <- NA_real_
        obj <- NA_real_
        deg <- ""
        lam <- ""

        rr <- tryCatch(
          run_kr(d, cfg),
          error = function(e) {
            status <<- "error"
            err <<- conditionMessage(e)
            NULL
          }
        )
        if (!is.null(rr)) {
          el <- rr$elapsed
          obj <- rr$objective
          deg <- rr$degree
          lam <- rr$lambda
        }

        rows[[idx]] <- data.frame(
          scenario = scn,
          seed_policy = spol,
          seed = as.integer(sd),
          profile = pname,
          nmulti = as.integer(cfg$nmulti),
          max_bb_eval = as.integer(cfg$max_bb_eval),
          status = status,
          elapsed_sec = el,
          objective = obj,
          degree = deg,
          lambda = lam,
          error = err,
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }
}

raw <- do.call(rbind, rows)
raw$replicate <- ave(
  seq_len(nrow(raw)),
  interaction(raw$scenario, raw$profile, raw$seed_policy, raw$seed, drop = TRUE),
  FUN = seq_along
)

ok <- raw[raw$status == "ok" & is.finite(raw$elapsed_sec) & is.finite(raw$objective), ]

cmp_rows <- list()
cidx <- 1L
tol <- 1e-6
for (scn in unique(ok$scenario)) {
  for (spol in unique(ok$seed_policy)) {
    base <- ok[
      ok$scenario == scn & ok$seed_policy == spol & ok$profile == "kr_current",
      c("seed", "replicate", "elapsed_sec", "objective")
    ]
    names(base)[3:4] <- c("elapsed_base", "objective_base")
    if (nrow(base) == 0) next

    for (pname in unique(ok$profile[ok$scenario == scn & ok$seed_policy == spol])) {
      cand <- ok[
        ok$scenario == scn & ok$seed_policy == spol & ok$profile == pname,
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
        worst_obj_diff = max(od, na.rm = TRUE),
        best_obj_diff = min(od, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      cidx <- cidx + 1L
    }
  }
}
cmp <- do.call(rbind, cmp_rows)

agg <- aggregate(cbind(worsen_count, improve_count, tie_count) ~ profile, data = cmp, FUN = sum)
cov <- aggregate(mean_elapsed_pct_change ~ profile, data = cmp, FUN = length)
names(cov)[2] <- "cells_present"
msp <- aggregate(mean_elapsed_pct_change ~ profile, data = cmp, FUN = mean)
names(msp)[2] <- "mean_elapsed_pct_change_overall"
mdp <- aggregate(median_elapsed_pct_change ~ profile, data = cmp, FUN = median)
names(mdp)[2] <- "median_elapsed_pct_change_overall"
woc <- aggregate(worst_obj_diff ~ profile, data = cmp, FUN = max)
names(woc)[2] <- "worst_obj_diff_overall"
boc <- aggregate(best_obj_diff ~ profile, data = cmp, FUN = min)
names(boc)[2] <- "best_obj_diff_overall"

best_seed <- aggregate(objective ~ scenario + seed_policy + seed + profile, data = ok, FUN = min)
best_baseline <- best_seed[best_seed$profile == "kr_current", c("scenario", "seed_policy", "seed", "objective")]
names(best_baseline)[4] <- "objective_base"
best_merge <- merge(best_seed, best_baseline, by = c("scenario", "seed_policy", "seed"), all.x = TRUE)
best_merge$best_obj_diff <- best_merge$objective - best_merge$objective_base
best_tb <- aggregate(cbind(best_better = as.integer(best_merge$best_obj_diff < -tol),
                           best_worse = as.integer(best_merge$best_obj_diff > tol)) ~ profile,
                     data = best_merge, FUN = sum)

out <- Reduce(function(x, y) merge(x, y, by = "profile", all = TRUE),
              list(agg, cov, msp, mdp, woc, boc, best_tb))

strict_safe <- out[
  out$cells_present == 4 & out$worsen_count == 0,
  ,
  drop = FALSE
]
strict_safe <- strict_safe[order(strict_safe$mean_elapsed_pct_change_overall), , drop = FALSE]

raw_path <- paste0(out_prefix, "_raw.csv")
cmp_path <- paste0(out_prefix, "_compare.csv")
agg_path <- paste0(out_prefix, "_agg.csv")
safe_path <- paste0(out_prefix, "_strict_safe.csv")

write.csv(raw, raw_path, row.names = FALSE)
write.csv(cmp, cmp_path, row.names = FALSE)
write.csv(out, agg_path, row.names = FALSE)
write.csv(strict_safe, safe_path, row.names = FALSE)

cat("WROTE_RAW:", raw_path, "\n")
cat("WROTE_COMPARE:", cmp_path, "\n")
cat("WROTE_AGG:", agg_path, "\n")
cat("WROTE_STRICT_SAFE:", safe_path, "\n")
