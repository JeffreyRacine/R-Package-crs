#!/usr/bin/env Rscript

suppressWarnings(suppressPackageStartupMessages(library(crs)))

make_data <- function(n, seed) {
  set.seed(seed)
  x <- runif(n)
  z <- factor(sample(c("A", "B"), n, replace = TRUE))
  y <- sin(2 * pi * x) + as.numeric(z == "B") + rnorm(n, sd = 0.2)
  list(
    x = x,
    z = z,
    y = y,
    xz = data.frame(x = x, z = z),
    df = data.frame(y = y, x = x, z = z)
  )
}

safe_collapse <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return("")
  }
  paste(as.character(x), collapse = ";")
}

run_fr <- function(dat, cfg) {
  out <- frscvNOMAD(
    xz = dat$xz,
    y = dat$y,
    complexity = "degree",
    segments = 2,
    degree.max = cfg$degree_max,
    degree.min = cfg$degree_min,
    nmulti = cfg$nmulti,
    max.bb.eval = cfg$max_bb_eval,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )
  list(
    objective = unname(out$cv.objc[1]),
    degree = safe_collapse(out$degree),
    segments = safe_collapse(out$segments),
    lambda = safe_collapse(out$lambda),
    bws = ""
  )
}

run_kr <- function(dat, cfg) {
  out <- krscvNOMAD(
    xz = dat$xz,
    y = dat$y,
    complexity = "degree",
    segments = 2,
    degree.max = cfg$degree_max,
    degree.min = cfg$degree_min,
    nmulti = cfg$nmulti,
    max.bb.eval = cfg$max_bb_eval,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )
  list(
    objective = unname(out$cv.objc[1]),
    degree = safe_collapse(out$degree),
    segments = safe_collapse(out$segments),
    lambda = safe_collapse(out$lambda),
    bws = ""
  )
}

run_glp <- function(dat, cfg) {
  out <- npglpreg(
    y ~ x + z,
    data = dat$df,
    cv = "degree-bandwidth",
    nmulti = cfg$nmulti,
    max.bb.eval = cfg$max_bb_eval_glp,
    degree.max = cfg$degree_max,
    degree.min = cfg$degree_min,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )
  list(
    objective = unname(out$fv[1]),
    degree = safe_collapse(out$degree),
    segments = "",
    lambda = "",
    bws = safe_collapse(signif(out$bws, digits = 12))
  )
}

run_snomadr_basic <- function(cfg) {
  eval.f <- function(x) {
    f <- c(Inf, Inf, Inf)
    n <- length(x)

    if (n == 5 && (is.double(x) || is.integer(x))) {
      f[1] <- x[5]
      f[2] <- sum((x - 1)^2) - 25
      f[3] <- 25 - sum((x + 1)^2)
    }

    as.double(f)
  }

  x0 <- rep(0.0, 5)
  bbin <- c(1, 1, 1, 1, 1)
  lb <- rep(-6.0, 5)
  ub <- c(5.0, 6.0, 7.0, 1000000, 100000)
  bbout <- c(0, 2, 1)
  opts <- list(
    "MAX_BB_EVAL" = as.integer(cfg$max_bb_eval_basic),
    "MIN_MESH_SIZE" = 0.001,
    "INITIAL_MESH_SIZE" = 0.1,
    "MIN_FRAME_SIZE" = 1,
    "DISPLAY_DEGREE" = 0
  )

  out <- snomadr(
    eval.f = eval.f,
    n = 5,
    x0 = x0,
    bbin = bbin,
    bbout = bbout,
    lb = lb,
    ub = ub,
    opts = opts,
    display.nomad.progress = FALSE
  )

  list(
    objective = unname(out$objective[1]),
    degree = "",
    segments = "",
    lambda = "",
    bws = safe_collapse(signif(out$solution, digits = 12))
  )
}

run_case <- function(case_id, seed, seed_policy, cfg) {
  dat <- make_data(cfg$n, seed)
  warning_msgs <- character(0)
  status <- "ok"
  err_msg <- ""
  result <- list(
    objective = NA_real_,
    degree = "",
    segments = "",
    lambda = "",
    bws = ""
  )
  elapsed <- as.numeric(system.time({
    result <- withCallingHandlers(
      tryCatch(
        switch(
          case_id,
          frscvNOMAD = run_fr(dat, cfg),
          krscvNOMAD = run_kr(dat, cfg),
          npglpreg = run_glp(dat, cfg),
          snomadr_basic_lib = run_snomadr_basic(cfg),
          stop(paste("unknown case", case_id))
        ),
        error = function(e) {
          status <<- "error"
          err_msg <<- conditionMessage(e)
          list(
            objective = NA_real_,
            degree = "",
            segments = "",
            lambda = "",
            bws = ""
          )
        }
      ),
      warning = function(w) {
        warning_msgs <<- c(warning_msgs, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
  })[["elapsed"]])

  param_signature <- paste(
    paste0("degree=", result$degree),
    paste0("segments=", result$segments),
    paste0("lambda=", result$lambda),
    paste0("bws=", result$bws),
    sep = "|"
  )

  data.frame(
    case = case_id,
    seed_policy = seed_policy,
    seed = as.integer(seed),
    n = as.integer(cfg$n),
    nmulti = as.integer(cfg$nmulti),
    max_bb_eval = as.integer(cfg$max_bb_eval),
    elapsed_sec = elapsed,
    objective = as.numeric(result$objective),
    degree = result$degree,
    segments = result$segments,
    lambda = result$lambda,
    bws = result$bws,
    param_signature = param_signature,
    status = status,
    warnings = safe_collapse(unique(warning_msgs)),
    error = err_msg,
    stringsAsFactors = FALSE
  )
}

summarize_raw <- function(raw_df) {
  ok <- raw_df[raw_df$status == "ok", , drop = FALSE]
  if (nrow(ok) == 0) {
    return(data.frame(
      case = character(0),
      seed_policy = character(0),
      runs = integer(0),
      mean_elapsed_sec = numeric(0),
      median_elapsed_sec = numeric(0),
      mean_objective = numeric(0),
      median_objective = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  key <- paste(ok$case, ok$seed_policy, sep = "||")
  parts <- split(ok, key)
  out <- lapply(parts, function(df) {
    data.frame(
      case = df$case[1],
      seed_policy = df$seed_policy[1],
      runs = nrow(df),
      mean_elapsed_sec = mean(df$elapsed_sec, na.rm = TRUE),
      median_elapsed_sec = stats::median(df$elapsed_sec, na.rm = TRUE),
      mean_objective = mean(df$objective, na.rm = TRUE),
      median_objective = stats::median(df$objective, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

parity_raw <- function(raw_df) {
  fixed <- raw_df[raw_df$status == "ok" & raw_df$seed_policy == "fixed_seed", , drop = FALSE]
  if (nrow(fixed) == 0) {
    return(data.frame(
      case = character(0),
      runs = integer(0),
      max_abs_objective_diff = numeric(0),
      unique_param_signatures = integer(0),
      stringsAsFactors = FALSE
    ))
  }
  parts <- split(fixed, fixed$case)
  out <- lapply(parts, function(df) {
    base_obj <- df$objective[1]
    data.frame(
      case = df$case[1],
      runs = nrow(df),
      max_abs_objective_diff = max(abs(df$objective - base_obj), na.rm = TRUE),
      unique_param_signatures = length(unique(df$param_signature)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

run_suite <- function(cfg, out_prefix) {
  if (is.null(cfg$max_bb_eval_basic)) {
    cfg$max_bb_eval_basic <- cfg$max_bb_eval
  }
  cases <- c("frscvNOMAD", "krscvNOMAD", "npglpreg", "snomadr_basic_lib")
  rows <- list()
  idx <- 1L
  for (seed in cfg$fixed_seeds) {
    for (case_id in cases) {
      rows[[idx]] <- run_case(case_id, seed, "fixed_seed", cfg)
      idx <- idx + 1L
    }
  }
  for (seed in cfg$varying_seeds) {
    for (case_id in cases) {
      rows[[idx]] <- run_case(case_id, seed, "varying_seed", cfg)
      idx <- idx + 1L
    }
  }
  raw_df <- do.call(rbind, rows)
  summary_df <- summarize_raw(raw_df)
  parity_df <- parity_raw(raw_df)

  raw_path <- paste0(out_prefix, "_raw.csv")
  summary_path <- paste0(out_prefix, "_summary.csv")
  parity_path <- paste0(out_prefix, "_parity.rds")

  write.csv(raw_df, raw_path, row.names = FALSE)
  write.csv(summary_df, summary_path, row.names = FALSE)
  saveRDS(parity_df, parity_path)

  cat("WROTE_RAW:", raw_path, "\n")
  cat("WROTE_SUMMARY:", summary_path, "\n")
  cat("WROTE_PARITY:", parity_path, "\n")
}
