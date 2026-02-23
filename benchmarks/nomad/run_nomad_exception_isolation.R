#!/usr/bin/env Rscript

suppressWarnings(suppressPackageStartupMessages(library(crs)))

args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(file_arg) > 0) normalizePath(sub("^--file=", "", file_arg[1])) else NA_character_

mk_data <- function(scn, n, seed) {
  set.seed(seed)
  if (scn == "simple") {
    x <- runif(n)
    z <- factor(sample(c("A", "B"), n, replace = TRUE))
    y <- sin(2 * pi * x) + as.numeric(z == "B") + rnorm(n, sd = 0.2)
    return(list(df = data.frame(y = y, x = x, z = z), formula = y ~ x + z))
  }
  x1 <- runif(n)
  x2 <- runif(n)
  z <- factor(sample(c("A", "B", "C"), n, replace = TRUE, prob = c(0.35, 0.35, 0.30)))
  y <- sin(2 * pi * x1) + cos(2 * pi * x2) + 0.5 * x1 * x2 +
    as.numeric(z == "B") * 0.4 - as.numeric(z == "C") * 0.25 + rnorm(n, sd = 0.2)
  list(df = data.frame(y = y, x1 = x1, x2 = x2, z = z), formula = y ~ x1 + x2 + z)
}

profile_opts <- function(profile) {
  switch(profile,
    baseline = list(),
    cs_opt = list("CS_OPTIMIZATION" = "yes"),
    nm_opt = list("NM_OPTIMIZATION" = "yes"),
    qp_opt = list("QP_OPTIMIZATION" = "yes"),
    random_opt = list("RANDOM_ALGO_OPTIMIZATION" = "yes"),
    noquad = list("QUAD_MODEL_SEARCH" = "no", "EVAL_QUEUE_SORT" = "DIR_LAST_SUCCESS"),
    quadbox1 = list("QUAD_MODEL_SEARCH_BOX_FACTOR" = "1.0", "QUAD_MODEL_BOX_FACTOR" = "1.0"),
    disco_opt = list("DISCO_MADS_OPTIMIZATION" = "yes"),
    vns_opt = list("VNS_MADS_OPTIMIZATION" = "yes"),
    stop(sprintf("unknown profile: %s", profile))
  )
}

run_npglpreg_once <- function(scenario, seed, n, profile) {
  d <- mk_data(scenario, n, seed)
  opts <- profile_opts(profile)
  status <- "ok"
  err <- ""
  obj <- NA_real_
  elapsed <- NA_real_
  tryCatch({
    elapsed <- as.numeric(system.time({
      fit <- npglpreg(
        formula = d$formula,
        data = d$df,
        cv = "degree-bandwidth",
        nmulti = 2,
        max.bb.eval = 60,
        degree.max = 4,
        degree.min = 0,
        display.warnings = FALSE,
        display.nomad.progress = FALSE,
        opts = opts
      )
    })[["elapsed"]])
    obj <- as.numeric(fit$fv[1])
  }, error = function(e) {
    status <<- "error"
    err <<- conditionMessage(e)
  })
  data.frame(
    scenario = scenario,
    seed = as.integer(seed),
    n = as.integer(n),
    profile = profile,
    status = status,
    elapsed_sec = elapsed,
    objective = obj,
    error = gsub("[\r\n]+", " | ", err),
    stringsAsFactors = FALSE
  )
}

# Single-shot mode for subprocess-isolated call.
if (length(args) >= 1 && identical(args[[1]], "--single")) {
  if (length(args) < 5) {
    stop("usage: run_nomad_exception_isolation.R --single <scenario> <seed> <n> <profile>")
  }
  out <- run_npglpreg_once(args[[2]], as.integer(args[[3]]), as.integer(args[[4]]), args[[5]])
  out$mode <- "isolated"
  out <- out[, c("mode", "scenario", "seed", "n", "profile", "status", "elapsed_sec", "objective", "error")]
  write.table(out, file = stdout(), sep = ",", row.names = FALSE, col.names = TRUE, quote = TRUE)
  quit(save = "no", status = 0)
}

out_prefix <- if (length(args) >= 1) args[[1]] else "/tmp/crs_nomad_exception_isolation"
n <- if (length(args) >= 2) as.integer(args[[2]]) else 180L

if (!is.finite(n) || n <= 0) {
  stop("n must be a positive integer")
}
if (is.na(script_path) || !nzchar(script_path) || !file.exists(script_path)) {
  stop("Unable to resolve script path for isolated subprocess mode")
}

run_isolated <- function(scenario, seed, n, profile) {
  cmd <- c(script_path, "--single", scenario, as.character(seed), as.character(n), profile)
  out <- tryCatch(system2("Rscript", cmd, stdout = TRUE, stderr = TRUE), error = function(e) character(0))
  if (length(out) == 0) {
    return(data.frame(
      mode = "isolated",
      scenario = scenario,
      seed = as.integer(seed),
      n = as.integer(n),
      profile = profile,
      status = "error",
      elapsed_sec = NA_real_,
      objective = NA_real_,
      error = "system2 failed",
      stringsAsFactors = FALSE
    ))
  }
  header_idx <- which(grepl("^\"mode\",\"scenario\"", out))[1]
  if (is.na(header_idx)) {
    return(data.frame(
      mode = "isolated",
      scenario = scenario,
      seed = as.integer(seed),
      n = as.integer(n),
      profile = profile,
      status = "error",
      elapsed_sec = NA_real_,
      objective = NA_real_,
      error = paste(out, collapse = " | "),
      stringsAsFactors = FALSE
    ))
  }
  txt <- paste(out[header_idx:length(out)], collapse = "\n")
  read.csv(text = txt, stringsAsFactors = FALSE)
}

scenarios <- c("simple", "hard")
seeds <- c(42L, 1001L)
sequence <- c("baseline", "disco_opt", "cs_opt", "nm_opt", "baseline", "vns_opt", "baseline")

rows <- list()
idx <- 1L
for (scn in scenarios) {
  for (sd in seeds) {
    for (p in sequence) {
      x <- run_npglpreg_once(scn, sd, n, p)
      x$mode <- "in_session"
      x <- x[, c("mode", "scenario", "seed", "n", "profile", "status", "elapsed_sec", "objective", "error")]
      rows[[idx]] <- x
      idx <- idx + 1L
    }
    for (p in sequence) {
      rows[[idx]] <- run_isolated(scn, sd, n, p)
      idx <- idx + 1L
    }
  }
}
raw <- do.call(rbind, rows)

summ_rows <- list()
j <- 1L
for (scn in scenarios) {
  for (sd in seeds) {
    sub_in <- raw[raw$mode == "in_session" & raw$scenario == scn & raw$seed == sd, ]
    sub_iso <- raw[raw$mode == "isolated" & raw$scenario == scn & raw$seed == sd, ]
    pre_in <- sub_in$objective[which(sub_in$profile == "baseline")[1]]
    post_in <- sub_in$objective[which(sub_in$profile == "baseline")[2]]
    post2_in <- sub_in$objective[which(sub_in$profile == "baseline")[3]]
    pre_iso <- sub_iso$objective[which(sub_iso$profile == "baseline")[1]]
    post_iso <- sub_iso$objective[which(sub_iso$profile == "baseline")[2]]
    post2_iso <- sub_iso$objective[which(sub_iso$profile == "baseline")[3]]

    summ_rows[[j]] <- data.frame(
      scenario = scn,
      seed = as.integer(sd),
      baseline_delta_after_disco_in = post_in - pre_in,
      baseline_delta_after_disco_isolated = post_iso - pre_iso,
      baseline_delta_after_vns_in = post2_in - pre_in,
      baseline_delta_after_vns_isolated = post2_iso - pre_iso,
      cs_status_in = sub_in$status[sub_in$profile == "cs_opt"][1],
      cs_status_isolated = sub_iso$status[sub_iso$profile == "cs_opt"][1],
      nm_status_in = sub_in$status[sub_in$profile == "nm_opt"][1],
      nm_status_isolated = sub_iso$status[sub_iso$profile == "nm_opt"][1],
      stringsAsFactors = FALSE
    )
    j <- j + 1L
  }
}
summary_df <- do.call(rbind, summ_rows)

raw_path <- paste0(out_prefix, "_raw.csv")
summary_path <- paste0(out_prefix, "_summary.csv")
write.csv(raw, raw_path, row.names = FALSE)
write.csv(summary_df, summary_path, row.names = FALSE)

cat("WROTE_RAW:", raw_path, "\n")
cat("WROTE_SUMMARY:", summary_path, "\n")
print(summary_df)
