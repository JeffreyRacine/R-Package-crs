args <- commandArgs(trailingOnly = TRUE)
label <- if (length(args) >= 1L) args[[1L]] else "candidate"
out_dir <- if (length(args) >= 2L) args[[2L]] else file.path(tempdir(), "kernel_native_timing")
n <- if (length(args) >= 3L) as.integer(args[[3L]]) else 5000L
reps <- if (length(args) >= 4L) as.integer(args[[4L]]) else 100L

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(crs))

set.seed(20260611L)
x <- cbind(runif(n), runif(n))
z <- cbind(
  z1 = sample(1:12, n, replace = TRUE),
  z2 = sample(1:10, n, replace = TRUE),
  z3 = sample(1:6, n, replace = TRUE)
)
is_ordered <- c(FALSE, TRUE, FALSE)
lambda <- c(0.35, 0.72, 0.45)
z_unique <- crs:::uniquecombs(z)
ind <- attr(z_unique, "index")
ind_vals <- unique(ind)
ind_list <- lapply(ind_vals, function(v) ind == v)
nrow_z_unique <- NROW(z_unique)
weights <- 0.4 + runif(n)

signal <- sin(2 * pi * x[, 1L]) + 0.7 * cos(2 * pi * x[, 2L]) +
  z[, 1L] / max(z[, 1L]) + 0.35 * z[, 2L] / max(z[, 2L]) -
  0.25 * z[, 3L] / max(z[, 3L])
y <- signal + rnorm(n, sd = 0.25 * sd(signal))
K <- matrix(c(4, 2, 3, 2), ncol = 2L, byrow = TRUE)

time_loop <- function(expr_fun, reps) {
  gc(FALSE)
  invisible(expr_fun())
  elapsed <- system.time({
    for (i in seq_len(reps)) invisible(expr_fun())
  })[["elapsed"]]
  elapsed / reps
}

kernel_sweep <- function() {
  out <- vector("list", length(ind_vals))
  for (i in seq_along(ind_vals)) {
    out[[i]] <- crs:::prod.kernel.matrix(
      Z = z,
      z = z_unique[ind_vals[i], ],
      lambda = lambda,
      is.ordered.z = is_ordered
    )
  }
  out
}

run_cv <- function(basis, weighted, cv.func) {
  crs:::cv.kernel.spline(
    x = x,
    y = y,
    z = z,
    K = K,
    lambda = lambda,
    z.unique = z_unique,
    ind = ind,
    ind.vals = ind_vals,
    ind.list = ind_list,
    nrow.z.unique = nrow_z_unique,
    is.ordered.z = is_ordered,
    knots = "quantiles",
    basis = basis,
    cv.func = cv.func,
    weights = if (weighted) weights else NULL,
    tau = NULL,
    singular.ok = FALSE,
    display.warnings = FALSE,
    use.ridge = FALSE,
    smooth.penalty = TRUE,
    penalty.scale = 1000,
    use.gram.cv = TRUE,
    gram.rcond.min = 1e-8,
    record.gram.stats = FALSE,
    use.cell.cache = TRUE
  )
}

grid <- expand.grid(
  basis = c("additive", "tensor", "glp"),
  weighted = c(FALSE, TRUE),
  cv.func = c("cv.ls", "cv.gcv", "cv.aic"),
  stringsAsFactors = FALSE
)

rows <- vector("list", nrow(grid) + 1L)
rows[[1L]] <- data.frame(
  label = label,
  n = n,
  reps = reps,
  basis = "kernel_sweep",
  weighted = NA,
  cv.func = NA,
  objective = NA_real_,
  per_eval = time_loop(kernel_sweep, reps),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(grid))) {
  basis <- grid$basis[i]
  weighted <- grid$weighted[i]
  cv.func <- grid$cv.func[i]
  objective <- run_cv(basis, weighted, cv.func)
  rows[[i + 1L]] <- data.frame(
    label = label,
    n = n,
    reps = reps,
    basis = basis,
    weighted = weighted,
    cv.func = cv.func,
    objective = as.numeric(objective),
    per_eval = time_loop(function() run_cv(basis, weighted, cv.func), reps),
    stringsAsFactors = FALSE
  )
}

results <- do.call(rbind, rows)
csv_path <- file.path(out_dir, paste0("kernel_native_fixed_eval_", label, ".csv"))
write.csv(results, csv_path, row.names = FALSE)
cat("wrote", csv_path, "\n")
