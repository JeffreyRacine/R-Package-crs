# Codex Work Summary (crs package)

## Archive Note (2026-03-13)

This summary is archived as historical implementation provenance.

Reason:

1. it records past investigations and keep/reject decisions,
2. it is not the current operational guide for `crs`,
3. retaining it in archive keeps the repo root focused on live package docs.

This document summarizes the investigations and changes considered since Codex installation, including what was kept vs. not kept and why. It also records benchmark highlights where available.

## 1) Native hat diagonal + .lm.fit path
**Idea:** Replace `lm.fit + hatvalues()` with `.lm.fit + hat.from.lm.fit()` and add native leverage computation to reduce R overhead.

**Changes kept:**
- Added native leverage routine `crs_hat_diag()` (C) to compute hat diagonal using `dqrqy` in C.
- Registered in `crs_init.c` and used in `hat.from.lm.fit` (via `spline.R`).
- `Makevars` updated to include `hat_diag.o`.

**Why kept:**
- Microbenchmarks and end-to-end CV runs showed consistent improvements at moderate/large sizes.
- Combined pipeline `.lm.fit + hat.from.lm.fit` is faster than `lm.fit + hat(fit$qr)` and `lm + hatvalues()`.

**Key timings (n=10000, M=100, k=10/50/100):**
- `.lm.fit + hat.from.lm.fit` vs `lm + hatvalues` (full pipeline):
  - k=10: ~2.05x faster (mean)
  - k=50: ~1.32x faster (mean)
  - k=100: ~1.14x faster (mean)
- `.lm.fit + hat.from.lm.fit` vs `lm.fit + hat(fit$qr)`:
  - k=10: ~1.08x faster (mean)
  - k=50: ~1.18x faster (mean)
  - k=100: ~1.16x faster (mean)

**Notes:**
- Precision differences are negligible in practice; output equivalence holds within expected numeric tolerance.

## 2) Cache group index vectors / avoid repeated as.matrix
**Idea:** Reduce repeated group/reshape operations (e.g., cache `ind.list` and avoid repeated `as.matrix`).

**Changes kept:**
- Step 1: Cache group index vectors (`ind.list`).
- Step 2: Avoid repeated `as.matrix` conversions in hot paths.

**Why kept:**
- Microbenchmarks showed small but consistent improvements (close to ~1% or more), enough to justify cleanup.

## 3) CV penalty: data-driven defaults
**Idea:** Replace `cv.maxPenalty = sqrt(.Machine$double.xmax)` with a data-driven penalty tied to the intercept-only CV baseline.

**Changes kept:**
- Added `resolve_cv_maxPenalty()` to compute a penalty from the intercept-only model.
- Updated call sites to pass `cv.func` so `cv.ls/cv.gcv` and `cv.aic` are scaled appropriately.
  - `cv.ls/cv.gcv`: `10 * LOO-CV(intercept-only)`
  - `cv.aic`: `AICc(intercept-only) + 10`
- Added unit tests in `tests/testthat/test-cv-maxPenalty.R`.

**Why kept:**
- More stable and interpretable than machine-max penalties.
- Better aligned with objective scale (log/MSE for AICc vs. LOO MSE for CV).

**Notes:**
- Weighted versions are handled using weighted mean/variance.

## 4) `resolve_cv_maxPenalty()` placement + utility helpers
**Idea:** Centralize shared helpers in `R/util.R` and remove duplicate definitions in `np.regression.glp.R`.

**Changes kept:**
- Removed fallback `NZD` in `np.regression.glp.R` (util.R version is more robust).
- Moved `is.fullrank()` and `scale_robust()` into `util.R`.

**Why kept:**
- Ensures a single authoritative implementation and avoids ambiguity.

## 5) `NZD_pos` guard for nonnegative denominators
**Idea:** Use `NZD_pos` when kernel sums are guaranteed nonnegative (`ckerorder == 2`), otherwise use `NZD`.

**Changes kept:**
- Added `NZD_den <- if (ckerorder == 2) NZD_pos else NZD` in:
  - `glpregEst()`
  - `minimand.cv.ls()`
  - `minimand.cv.aic()`
- Replaced `NZD(...)` with `NZD_den(...)` in GLP/LOO/AIC paths.

**Why kept:**
- Empirical check: denominators are positive for `ckerorder=2` but can be negative for `ckerorder>2`.
- Guard avoids bias for higher-order kernels while allowing the faster/simpler positive-only handling for default order.

**Empirical check (degree>0, synthetic data):**
- `ckerorder=2`: denominators positive.
- `ckerorder=4/6`: denominators can be negative.

## 6) Items investigated but NOT kept
**A) Caching inside `cv.kernel.spline`**
- Tried caching within `cv.kernel.spline`/`prod.spline` (basis reuse).
- **Rejected**: measured slowdown vs. baseline in microbenchmarks.

**B) Candidate 4: `npglpreg` CV loop caching**
- Tried caching `scale_robust()`/numeric indices and passing into minimands.
- **Rejected**: no speedup; slightly slower in the synthetic check.

## 7) Pre/post package timing comparisons
**Fixed-parameter comparisons (degree/segments/bandwidth fixed)**
- Ensured output equivalence (max abs diff = 0).
- Large gains in `crs`, `crsiv`, `crsivderiv`, modest slowdown in `npglpreg` with fixed bandwidths.

**Default comparisons**
- Default CV searches can differ in selected models due to mixed-integer optimization and different evaluation paths.
- Differences are expected and not interpreted as regressions without fixed-parameter runs.

## 8) Tests and checks
- Added unit tests for `resolve_cv_maxPenalty` (cv.ls/cv.gcv/cv.aic; weighted/unweighted).
- `R CMD build` and `R CMD check` pass (network warnings only).

---

## Notes on reproducibility
- CV/search is stochastic and sensitive to floating-point differences; fixed parameters are required for strict equivalence.
- When evaluating performance, use fixed tuning or large enough runs to reduce timing noise.
