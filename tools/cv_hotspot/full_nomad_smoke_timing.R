args <- commandArgs(trailingOnly = TRUE)
label <- if (length(args) >= 1L) args[[1L]] else "candidate"
out_dir <- if (length(args) >= 2L) args[[2L]] else file.path(tempdir(), "full_nomad_smoke")
n <- if (length(args) >= 3L) as.integer(args[[3L]]) else 900L
nseed <- if (length(args) >= 4L) as.integer(args[[4L]]) else 5L

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(crs))
options(crs.messages = FALSE)

make_data <- function(n) {
  set.seed(20260611L)
  x1 <- runif(n)
  x2 <- runif(n)
  z1 <- factor(sample(letters[1:4], n, replace = TRUE))
  z2 <- ordered(sample(1:5, n, replace = TRUE))
  signal <- sin(2 * pi * x1) + 0.65 * cos(2 * pi * x2) +
    as.numeric(z1) / max(as.numeric(z1)) -
    0.35 * as.numeric(z2) / max(as.numeric(z2))
  data.frame(
    y = signal + rnorm(n, sd = 0.35 * sd(signal)),
    x1 = x1,
    x2 = x2,
    z1 = z1,
    z2 = z2
  )
}

extract_scalar <- function(x, names) {
  for (nm in names) {
    if (!is.null(x[[nm]]) && length(x[[nm]]) > 0L) return(as.numeric(x[[nm]][1L]))
  }
  NA_real_
}

dat <- make_data(n)
seeds <- seq_len(nseed) + 88000L
rows <- vector("list", length(seeds))

for (i in seq_along(seeds)) {
  seed <- seeds[i]
  set.seed(seed)
  elapsed <- system.time({
    fit <- crs(
      y ~ x1 + x2 + z1 + z2,
      data = dat,
      kernel = TRUE,
      cv = "nomad",
      cv.threshold = 0,
      nmulti = 1,
      random.seed = seed,
      opts = list("MAX_BB_EVAL" = 160, "DISPLAY_DEGREE" = 0)
    )
  })[["elapsed"]]

  rows[[i]] <- data.frame(
    label = label,
    n = n,
    seed = seed,
    elapsed = elapsed,
    cv.score = extract_scalar(fit, c("cv.score", "cv")),
    degree = paste(fit$degree, collapse = "/"),
    segments = paste(fit$segments, collapse = "/"),
    lambda = paste(signif(fit$lambda, 8), collapse = "/"),
    basis = as.character(fit$basis)[1L],
    stringsAsFactors = FALSE
  )
}

results <- do.call(rbind, rows)
csv_path <- file.path(out_dir, paste0("full_nomad_smoke_", label, ".csv"))
write.csv(results, csv_path, row.names = FALSE)
cat("wrote", csv_path, "\n")
print(results)
