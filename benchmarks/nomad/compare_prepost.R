#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
pre_raw <- if (length(args) >= 1) args[[1]] else "/tmp/crs_nomad_pre_bench_raw.csv"
post_raw <- if (length(args) >= 2) args[[2]] else "/tmp/crs_nomad_post_bench_raw.csv"
perf_out <- if (length(args) >= 3) args[[3]] else "/tmp/crs_nomad_prepost_perf_summary.csv"
parity_out <- if (length(args) >= 4) args[[4]] else "/tmp/crs_nomad_prepost_parity_summary.csv"

pre <- read.csv(pre_raw, stringsAsFactors = FALSE)
post <- read.csv(post_raw, stringsAsFactors = FALSE)

join_keys <- c("case", "seed_policy", "seed")
keep_cols <- c(join_keys, "elapsed_sec", "objective", "degree", "segments", "lambda", "bws", "status")

pre <- pre[keep_cols]
post <- post[keep_cols]
colnames(pre)[!(colnames(pre) %in% join_keys)] <- paste0(colnames(pre)[!(colnames(pre) %in% join_keys)], "_pre")
colnames(post)[!(colnames(post) %in% join_keys)] <- paste0(colnames(post)[!(colnames(post) %in% join_keys)], "_post")

merged <- merge(pre, post, by = join_keys, all = FALSE, sort = FALSE)
merged$elapsed_pct_change <- 100 * (merged$elapsed_sec_post - merged$elapsed_sec_pre) / pmax(merged$elapsed_sec_pre, .Machine$double.eps)
merged$objective_abs_diff <- abs(merged$objective_post - merged$objective_pre)
merged$params_changed <- with(merged, degree_pre != degree_post | segments_pre != segments_post | lambda_pre != lambda_post | bws_pre != bws_post)

ok <- merged[merged$status_pre == "ok" & merged$status_post == "ok", , drop = FALSE]

key <- paste(ok$case, ok$seed_policy, sep = "||")
parts <- split(ok, key)

perf_rows <- lapply(parts, function(df) {
  data.frame(
    case = df$case[1],
    seed_policy = df$seed_policy[1],
    runs = nrow(df),
    mean_elapsed_pct_change = mean(df$elapsed_pct_change, na.rm = TRUE),
    median_elapsed_pct_change = stats::median(df$elapsed_pct_change, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
})

parity_rows <- lapply(parts, function(df) {
  data.frame(
    case = df$case[1],
    seed_policy = df$seed_policy[1],
    runs = nrow(df),
    max_abs_objective_diff = max(df$objective_abs_diff, na.rm = TRUE),
    params_changed_runs = sum(df$params_changed, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
})

perf_df <- do.call(rbind, perf_rows)
parity_df <- do.call(rbind, parity_rows)

write.csv(perf_df, perf_out, row.names = FALSE)
write.csv(parity_df, parity_out, row.names = FALSE)

cat("WROTE_PERF:", perf_out, "\n")
cat("WROTE_PARITY:", parity_out, "\n")
