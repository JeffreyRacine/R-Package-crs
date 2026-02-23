#!/usr/bin/env Rscript

suppressWarnings(suppressPackageStartupMessages(library(crs)))

args <- commandArgs(trailingOnly = TRUE)
out_prefix <- if (length(args) >= 1) args[[1]] else "/tmp/crs_nomad_dimmix_tradeoff"
n <- if (length(args) >= 2) as.integer(args[[2]]) else 100L

safe_collapse <- function(x) {
  if (is.null(x) || length(x) == 0) return("")
  paste(as.character(x), collapse = ";")
}

mk_data <- function(p_cont, p_cat, n, seed) {
  set.seed(seed)

  x <- as.data.frame(replicate(p_cont, runif(n), simplify = FALSE))
  names(x) <- paste0("x", seq_len(p_cont))

  z <- as.data.frame(lapply(seq_len(p_cat), function(j) {
    factor(
      sample(c("A", "B", "C"), n, replace = TRUE, prob = c(0.35, 0.35, 0.30)),
      levels = c("A", "B", "C")
    )
  }))
  names(z) <- paste0("z", seq_len(p_cat))

  y <- rep(0, n)
  if (p_cont >= 1) y <- y + sin(2 * pi * x[[1]])
  if (p_cont >= 2) y <- y + cos(2 * pi * x[[2]])
  if (p_cont >= 3) y <- y + 0.5 * x[[1]] * x[[3]]
  if (p_cont > 3) {
    for (j in 4:p_cont) y <- y + 0.1 * x[[j]]
  }

  for (j in seq_len(p_cat)) {
    zc <- z[[j]]
    y <- y + 0.25 * as.numeric(zc == "B") - 0.18 * as.numeric(zc == "C")
  }

  y <- y + rnorm(n, sd = 0.2)

  xz <- cbind(x, z)
  np_df <- cbind(data.frame(y = y), x, z)
  rhs <- paste(c(names(x), names(z)), collapse = " + ")
  np_formula <- as.formula(paste("y ~", rhs))

  list(
    xz = xz,
    y = y,
    np_df = np_df,
    np_formula = np_formula,
    segments = rep(2L, p_cont),
    degree_max = 4L,
    degree_min = 0L
  )
}

run_case <- function(case_id, dat, profile) {
  tm <- system.time({
    out <- switch(
      case_id,
      frscvNOMAD = frscvNOMAD(
        xz = dat$xz,
        y = dat$y,
        complexity = "degree",
        segments = dat$segments,
        degree.max = dat$degree_max,
        degree.min = dat$degree_min,
        nmulti = as.integer(profile$nmulti),
        max.bb.eval = as.integer(profile$max_bb_eval),
        opts = profile$opts,
        display.warnings = FALSE,
        display.nomad.progress = FALSE
      ),
      krscvNOMAD = krscvNOMAD(
        xz = dat$xz,
        y = dat$y,
        complexity = "degree",
        segments = dat$segments,
        degree.max = dat$degree_max,
        degree.min = dat$degree_min,
        nmulti = as.integer(profile$nmulti),
        max.bb.eval = as.integer(profile$max_bb_eval),
        opts = profile$opts,
        display.warnings = FALSE,
        display.nomad.progress = FALSE
      ),
      npglpreg = npglpreg(
        formula = dat$np_formula,
        data = dat$np_df,
        cv = "degree-bandwidth",
        nmulti = as.integer(profile$nmulti),
        max.bb.eval = as.integer(profile$max_bb_eval),
        degree.max = dat$degree_max,
        degree.min = dat$degree_min,
        opts = profile$opts,
        display.warnings = FALSE,
        display.nomad.progress = FALSE
      )
    )
  })[["elapsed"]]

  switch(
    case_id,
    frscvNOMAD = list(
      elapsed = as.numeric(tm),
      objective = out$cv.objc[1],
      signature = paste(
        paste0("deg=", safe_collapse(out$degree)),
        paste0("seg=", safe_collapse(out$segments)),
        paste0("inc=", safe_collapse(out$include)),
        sep = "|"
      )
    ),
    krscvNOMAD = list(
      elapsed = as.numeric(tm),
      objective = out$cv.objc[1],
      signature = paste(
        paste0("deg=", safe_collapse(out$degree)),
        paste0("seg=", safe_collapse(out$segments)),
        paste0("lam=", safe_collapse(signif(out$lambda, 10))),
        sep = "|"
      )
    ),
    npglpreg = list(
      elapsed = as.numeric(tm),
      objective = out$fv[1],
      signature = paste(
        paste0("deg=", safe_collapse(out$degree)),
        paste0("bws=", safe_collapse(signif(out$bws, 10))),
        sep = "|"
      )
    )
  )
}

profiles <- list(
  frscvNOMAD = list(
    fr_current = list(nmulti = 2L, max_bb_eval = 100L, opts = list()),
    fr_universal_aggr = list(
      nmulti = 2L,
      max_bb_eval = 140L,
      opts = list("DIRECTION_TYPE" = "ORTHO 2N")
    )
  ),
  krscvNOMAD = list(
    kr_current = list(nmulti = 2L, max_bb_eval = 100L, opts = list()),
    kr_nmtrial40 = list(
      nmulti = 2L,
      max_bb_eval = 100L,
      opts = list("NM_SEARCH_MAX_TRIAL_PTS_NFACTOR" = 40L)
    ),
    kr_universal_aggr = list(
      nmulti = 2L,
      max_bb_eval = 140L,
      opts = list("DIRECTION_TYPE" = "ORTHO 2N")
    )
  ),
  npglpreg = list(
    np_current = list(nmulti = 2L, max_bb_eval = 100L, opts = list()),
    np_universal_aggr = list(
      nmulti = 2L,
      max_bb_eval = 140L,
      opts = list("DIRECTION_TYPE" = "ORTHO 2N")
    ),
    np_noquad = list(
      nmulti = 2L,
      max_bb_eval = 100L,
      opts = list(
        "QUAD_MODEL_SEARCH" = "no",
        "EVAL_QUEUE_SORT" = "DIR_LAST_SUCCESS"
      )
    )
  )
)

baselines <- c(frscvNOMAD = "fr_current", krscvNOMAD = "kr_current", npglpreg = "np_current")
tols <- c(frscvNOMAD = 1e-10, krscvNOMAD = 1e-6, npglpreg = 5e-4)

dim_grid <- list(
  d1c1 = c(1L, 1L),
  d2c2 = c(2L, 2L),
  d3c3 = c(3L, 3L)
)

rows <- list()
idx <- 1L

for (dname in names(dim_grid)) {
  p_cont <- dim_grid[[dname]][1]
  p_cat <- dim_grid[[dname]][2]
  cat("DIM", dname, "p_cont", p_cont, "p_cat", p_cat, "\n")

  for (spol in c("fixed_seed", "varying_seed")) {
    seeds <- if (spol == "fixed_seed") rep(42L, 5L) else c(1001L, 1002L, 1003L, 1004L, 1005L)
    for (sd in seeds) {
      dat <- mk_data(p_cont = p_cont, p_cat = p_cat, n = n, seed = sd)
      for (case_id in names(profiles)) {
        for (pname in names(profiles[[case_id]])) {
          prof <- profiles[[case_id]][[pname]]
          status <- "ok"
          err <- ""
          elapsed <- NA_real_
          objective <- NA_real_
          sig <- ""

          rr <- tryCatch(
            run_case(case_id, dat, prof),
            error = function(e) {
              status <<- "error"
              err <<- conditionMessage(e)
              NULL
            }
          )
          if (!is.null(rr)) {
            elapsed <- rr$elapsed
            objective <- rr$objective
            sig <- rr$signature
          }

          rows[[idx]] <- data.frame(
            dim_id = dname,
            p_cont = as.integer(p_cont),
            p_cat = as.integer(p_cat),
            seed_policy = spol,
            seed = as.integer(sd),
            case = case_id,
            profile = pname,
            nmulti = as.integer(prof$nmulti),
            max_bb_eval = as.integer(prof$max_bb_eval),
            status = status,
            elapsed_sec = as.numeric(elapsed),
            objective = as.numeric(objective),
            signature = sig,
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
  interaction(raw$dim_id, raw$case, raw$profile, raw$seed_policy, raw$seed, drop = TRUE),
  FUN = seq_along
)

ok <- raw[raw$status == "ok" & is.finite(raw$elapsed_sec) & is.finite(raw$objective), ]

cmp_rows <- list()
cidx <- 1L
for (dname in unique(ok$dim_id)) {
  for (case_id in unique(ok$case)) {
    base_name <- baselines[[case_id]]
    tol <- tols[[case_id]]
    for (spol in unique(ok$seed_policy)) {
      base <- ok[
        ok$dim_id == dname & ok$case == case_id & ok$seed_policy == spol & ok$profile == base_name,
        c("seed", "replicate", "elapsed_sec", "objective")
      ]
      names(base)[3:4] <- c("elapsed_base", "objective_base")
      if (nrow(base) == 0) next

      for (pname in unique(ok$profile[ok$dim_id == dname & ok$case == case_id])) {
        cand <- ok[
          ok$dim_id == dname & ok$case == case_id & ok$seed_policy == spol & ok$profile == pname,
          c("seed", "replicate", "elapsed_sec", "objective")
        ]
        names(cand)[3:4] <- c("elapsed_cand", "objective_cand")
        if (nrow(cand) == 0) next

        m <- merge(base, cand, by = c("seed", "replicate"), all = FALSE)
        if (nrow(m) == 0) next

        dp <- 100 * (m$elapsed_cand - m$elapsed_base) / pmax(m$elapsed_base, .Machine$double.eps)
        od <- m$objective_cand - m$objective_base

        cmp_rows[[cidx]] <- data.frame(
          dim_id = dname,
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
          worst_obj_diff = max(od, na.rm = TRUE),
          best_obj_diff = min(od, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
        cidx <- cidx + 1L
      }
    }
  }
}
cmp <- do.call(rbind, cmp_rows)

agg_by_dim <- aggregate(
  cbind(worsen_count, improve_count, tie_count) ~ dim_id + case + profile,
  data = cmp,
  FUN = sum
)
cov_by_dim <- aggregate(mean_elapsed_pct_change ~ dim_id + case + profile, data = cmp, FUN = length)
names(cov_by_dim)[4] <- "cells_present"
mme_by_dim <- aggregate(mean_elapsed_pct_change ~ dim_id + case + profile, data = cmp, FUN = mean)
names(mme_by_dim)[4] <- "mean_elapsed_pct_change_overall"
med_by_dim <- aggregate(median_elapsed_pct_change ~ dim_id + case + profile, data = cmp, FUN = median)
names(med_by_dim)[4] <- "median_elapsed_pct_change_overall"
wod_by_dim <- aggregate(worst_obj_diff ~ dim_id + case + profile, data = cmp, FUN = max)
names(wod_by_dim)[4] <- "worst_obj_diff_overall"
bod_by_dim <- aggregate(best_obj_diff ~ dim_id + case + profile, data = cmp, FUN = min)
names(bod_by_dim)[4] <- "best_obj_diff_overall"

agg_dim <- Reduce(
  function(x, y) merge(x, y, by = c("dim_id", "case", "profile"), all = TRUE),
  list(agg_by_dim, cov_by_dim, mme_by_dim, med_by_dim, wod_by_dim, bod_by_dim)
)

agg_all <- aggregate(
  cbind(worsen_count, improve_count, tie_count) ~ case + profile,
  data = cmp,
  FUN = sum
)
cov_all <- aggregate(mean_elapsed_pct_change ~ case + profile, data = cmp, FUN = length)
names(cov_all)[3] <- "cells_present"
mme_all <- aggregate(mean_elapsed_pct_change ~ case + profile, data = cmp, FUN = mean)
names(mme_all)[3] <- "mean_elapsed_pct_change_overall"
med_all <- aggregate(median_elapsed_pct_change ~ case + profile, data = cmp, FUN = median)
names(med_all)[3] <- "median_elapsed_pct_change_overall"
wod_all <- aggregate(worst_obj_diff ~ case + profile, data = cmp, FUN = max)
names(wod_all)[3] <- "worst_obj_diff_overall"
bod_all <- aggregate(best_obj_diff ~ case + profile, data = cmp, FUN = min)
names(bod_all)[3] <- "best_obj_diff_overall"

agg_overall <- Reduce(
  function(x, y) merge(x, y, by = c("case", "profile"), all = TRUE),
  list(agg_all, cov_all, mme_all, med_all, wod_all, bod_all)
)

raw_path <- paste0(out_prefix, "_raw.csv")
cmp_path <- paste0(out_prefix, "_compare.csv")
dim_path <- paste0(out_prefix, "_agg_by_dim.csv")
ovr_path <- paste0(out_prefix, "_agg_overall.csv")

write.csv(raw, raw_path, row.names = FALSE)
write.csv(cmp, cmp_path, row.names = FALSE)
write.csv(agg_dim, dim_path, row.names = FALSE)
write.csv(agg_overall, ovr_path, row.names = FALSE)

cat("WROTE_RAW:", raw_path, "\n")
cat("WROTE_COMPARE:", cmp_path, "\n")
cat("WROTE_AGG_BY_DIM:", dim_path, "\n")
cat("WROTE_AGG_OVERALL:", ovr_path, "\n")
