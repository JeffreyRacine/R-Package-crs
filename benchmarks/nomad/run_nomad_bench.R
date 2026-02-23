#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
out_prefix <- if (length(args) >= 1) args[[1]] else "/tmp/crs_nomad_bench"

file_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
} else {
  script_dir <- normalizePath(getwd())
}

source(file.path(script_dir, "common_harness.R"))

cfg <- list(
  n = 120L,
  nmulti = 2L,
  degree_max = 4L,
  degree_min = 0L,
  max_bb_eval = 80L,
  max_bb_eval_glp = 80L,
  fixed_seeds = rep(42L, 5L),
  varying_seeds = c(1001L, 1002L, 1003L, 1004L, 1005L)
)

run_suite(cfg = cfg, out_prefix = out_prefix)
