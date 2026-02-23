#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
out_prefix <- if (length(args) >= 1) args[[1]] else "/tmp/crs_nomad_smoke"

file_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
} else {
  script_dir <- normalizePath(getwd())
}

source(file.path(script_dir, "common_harness.R"))

cfg <- list(
  n = 40L,
  nmulti = 1L,
  degree_max = 3L,
  degree_min = 0L,
  max_bb_eval = 30L,
  max_bb_eval_glp = 30L,
  fixed_seeds = rep(42L, 3L),
  varying_seeds = c(101L, 202L, 303L)
)

run_suite(cfg = cfg, out_prefix = out_prefix)
