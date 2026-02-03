## This function provides a matrix with all combinations of a vector K
## (containing degrees) and vector I (0/1 for exclude/include for
## factors).

# > print(mb)
# Unit: milliseconds
# expr      min       lq     mean   median       uq       max neval cld
# original 1.451195 1.479198 1.863492 1.501645 1.530489 50.018401   200   a
# optimized 1.393344 1.420445 1.749084 1.445844 1.502773  3.504967   200   a
#
# > fmt_report(mb)
#
# -------- Benchmark Summary --------
# Median:  1.04x  ( 96.3% as fast)  - optimized * is fastest by median.
# Mean:    1.07x  ( 93.9% as fast)  - optimized * is fastest by mean.

# ------------------------------------------------
# Optimized (identical)
# ------------------------------------------------
# Key changes:
#  - Build the argument list for expand.grid using vectorized list replication:
#      c(rep(list(K.vec1), num.x), rep(list(K.vec2), num.x), rep(list(0:1), num.z))
#    guarded exactly like the original to preserve column order and presence/absence.
#  - Tell expand.grid to skip computing out.attrs (small constant-factor savings):
#      KEEP.OUT.ATTRS = FALSE   (has no effect on the returned *matrix*).
#  - Everything else (including error check and column ordering) matches the original.
matrix.combn <- function(K.vec1, K.vec2 = NULL, num.x = 0, num.z = NULL) {
  if (num.x == 0 && num.z == 0) stop(" must provide at least one variable")

  # Construct the list of components in the exact same order as original:
  parts <- list()

  if (num.x > 0) {
    parts <- c(parts, rep(list(K.vec1), num.x))
    if (!is.null(K.vec2)) {
      parts <- c(parts, rep(list(K.vec2), num.x))
      if (!is.null(num.z)) parts <- c(parts, rep(list(0:1), num.z))
    } else {
      if (!is.null(num.z)) parts <- c(parts, rep(list(0:1), num.z))
    }
  } else {
    if (!is.null(num.z)) parts <- c(parts, rep(list(0:1), num.z))
  }

  if (length(parts) == 0L) {
    # Defensive: original would never reach here because of the stop() above.
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }

  as.matrix(do.call(expand.grid, c(parts, list(KEEP.OUT.ATTRS = FALSE))))
}
