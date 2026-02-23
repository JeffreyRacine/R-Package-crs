# CRS Modernization Plan (R + Non-NOMAD C)

Date: 2026-02-23  
Status: Planning only (no implementation started)

## Goal

Modernize `crs` to current best-practice R package engineering standards while preserving:

1. Functional behavior.
2. Numerical stability (except accepted optimizer-path drift).
3. Existing user-facing API/contracts.

Primary guidance basis:

1. `R Packages (2e)` (`r-pkgs`): package/release hygiene, dependency and check discipline.
2. `Advanced R (2e)` (`adv-r`): evaluation semantics, call construction, state management.
3. `R for Data Science (2e)` (`r4ds`): reproducible workflow, clear validation artifacts, communicable results.

## Hard Scope / Freeze Rules

Out of scope for this modernization project:

1. Embedded NOMAD source and C interface:
   - `/Users/jracine/Development/crs/src/nomad4_src/**`
   - `/Users/jracine/Development/crs/src/snomadr.cpp`
   - `/Users/jracine/Development/crs/src/snomadr.h`
2. NOMAD strategy/tuning/options behavior changes.
3. NOMAD docs/map files except cross-reference updates if needed.

In scope:

1. R-layer modernization (including NOMAD-adjacent wrappers only where changes are purely engineering and behavior-preserving).
2. Non-NOMAD C/C integration and safety hardening:
   - `RuniqueCombs.c`, `bspline.c`, `gsl_bspline.c`, `hat_diag.c`, `mgcv.c`, `crs_init.c`.

## What Was Learned From `np-master` / `np-npRmpi`

Modernization patterns that worked repeatedly and safely:

1. Small, tranche-based refactors (one risk class per commit).
2. Replace scalar footguns first:
   - `1:length(...)` -> `seq_along(...)`, `1:ncol(...)` -> `seq_len(ncol(...))`
   - scalar `ifelse(...)` -> scalar `if (...) ... else ...`
   - scalar `|`/`&` in controls -> `||`/`&&`
   - `options('x')$x` -> `getOption("x")` / `isTRUE(...)`
3. Reduce dynamic-eval/call-construction risk before deeper algorithmic rewrites.
4. Performance governance per checkpoint:
   - fixed-seed + varying-seed
   - mean/median deltas
   - parity summary
   - `/tmp` artifacts logged.
5. Tarball-first checks as authoritative gates.

## Baseline Audit Snapshot (`crs`)

Static inventory (current):

1. `.C(` callsites in `R/`: `3`
2. `.Call(` callsites in `R/`: `6`
3. `eval(parse(...))` in `R/`: `2` (both in `R/stepCV.R`)
4. `eval(` in `R/`: `7`
5. `do.call("<string>", ...)` in `R/`: `13`
6. `<<-` in `R/`: `15`
7. `ifelse(` in `R/`: `58`
8. `1:length(...)` in `R/`: `50`

Hotspot files (non-NOMAD-heavy):

1. `/Users/jracine/Development/crs/R/spline.R` (`ifelse` density).
2. `/Users/jracine/Development/crs/R/crs.R` (`eval(object$call$data)`, indexing patterns).
3. `/Users/jracine/Development/crs/R/crsiv.R` and `/Users/jracine/Development/crs/R/crsivderiv.R` (string `do.call`, repeated loop scaffolding).
4. `/Users/jracine/Development/crs/R/stepCV.R` (`eval(parse(...))`, duplicate model fitting in CV loop).
5. `/Users/jracine/Development/crs/R/console.R`, `/Users/jracine/Development/crs/R/frscv.R`, `/Users/jracine/Development/crs/R/krscv.R`, `/Users/jracine/Development/crs/R/np.regression.glp.R` (`<<-` usage).

Non-NOMAD C note identified during audit:

1. `/Users/jracine/Development/crs/src/gsl_bspline.c` has a likely leak path in `gsl_bspline_deriv` (`quantile_vec` allocation not freed; free call commented out).

## Prioritized Execution Plan (High ROI / Low Risk First)

### A/R1 - Low risk, high ROI (start here)

1. Scalar/index safety sweep in non-NOMAD R paths:
   - `1:length(...)` -> `seq_along(...)`
   - `1:ncol(...)` -> `seq_len(ncol(...))`
   - scalar `ifelse` -> scalar `if`.
2. Option-access hardening:
   - `options('crs.messages')$crs.messages` -> `isTRUE(getOption("crs.messages"))`.
3. String `do.call` cleanup:
   - `do.call("crs", ...)` -> `do.call(crs, ...)`
   - `do.call("expand.grid", ...)` -> `do.call(base::expand.grid, ...)`.
4. Super-assignment retirement where local state is sufficient:
   - `console`, `tmpPos`, `dots` mutation patterns.

Expected impact:

1. Lower regression risk from scalar edge cases.
2. Cleaner, faster control flow.
3. Better maintainability with minimal behavioral risk.

### A/R2 - Medium risk, high ROI

1. Remove `eval(parse(...))` in `R/stepCV.R` via structured formula construction (`reformulate`/`update.formula` path).
2. Remove redundant re-fit in `stepCV` CV scoring loop:
   - avoid `lm(y ~ X)` where equivalent quantities can be derived from existing fit/QR path.
3. Replace fragile `eval(object$call$data)` in `R/crs.R` with robust call/data resolution helper.
4. Reduce repeated call/list rebuilding in `crsiv*` loops by prebuilding immutable argument bundles.

Expected impact:

1. Performance gains in iterative/model-selection paths.
2. Reduced evaluation-semantics fragility.

### B/R2 - Medium risk, medium/high ROI

1. DRY refactor duplicated `crsiv`/`crsivderiv` loop scaffolding into shared internal helpers.
2. Hoist immutable computation out of iterative updates where possible:
   - precompute structures that do not vary across iterations (formula components, indexing maps, stable design metadata).
3. Normalize progress-console handling to side-effect-safe local state.

Expected impact:

1. Lower code duplication and maintenance risk.
2. Potential measurable runtime reductions in iterative routines.

### B/R1 - Low risk, medium ROI (non-NOMAD C hygiene)

1. Fix confirmed memory/safety issues in non-NOMAD C (`gsl_bspline.c` leak candidate).
2. Add/expand focused tests that exercise C wrappers (`gsl.bs`, `uniquecombs`, hat-diagonal path).
3. Verify registration and call signatures stay exact (`crs_init.c` gate).

### C/R3 - High risk, optional/deferred

1. `.C` -> `.Call` migration for legacy non-NOMAD entry points (`gsl_bspline`, `RuniqueCombs`) if and only if earlier tranches are stable.
2. Deep algorithmic rewrites inside large legacy routines (e.g., major `spline.R` internals) only after broad parity/performance signoff.

## Validation Gates (After Every Tranche)

Minimum required after each checkpoint:

1. Targeted `testthat` for touched surfaces.
2. Installed-package smoke scripts for impacted estimators.
3. Tarball-first:
   - from `/Users/jracine/Development`: `R CMD build crs`
   - from `/Users/jracine/Development`: `R CMD check --as-cran crs_<ver>.tar.gz`
4. For performance-affecting changes:
   - fixed-seed and varying-seed comparisons
   - mean/median deltas
   - parity summary (at least one objective/output parity check)
   - artifacts under `/tmp` with clear names.

Heavy example gate (periodic, not every commit):

1. `/Users/jracine/Development/crs/man/runcrs`
2. `R CMD build` + `R CMD check --as-cran`
3. `/Users/jracine/Development/crs/man/dontruncrs`
4. Compare `crs-Ex.Rout` and `crs-Ex.timings` versus baseline.

## Commit Strategy

1. One modernization theme per commit (small diff, narrow blast radius).
2. Commit message format: `modernize(<area>): <change>`.
3. Include validation artifact paths in commit body or companion note.
4. Local commits only (no push).

## Proposed First Implementation Slice (when authorized)

A/R1 tranche candidate:

1. `options('crs.messages')$...` -> `getOption`/`isTRUE` in:
   - `/Users/jracine/Development/crs/R/console.R`
   - `/Users/jracine/Development/crs/R/crs.R`
   - `/Users/jracine/Development/crs/R/clsd.R`
   - `/Users/jracine/Development/crs/R/np.regression.glp.R`
   - `/Users/jracine/Development/crs/R/zzz.R`
2. Retire low-risk string `do.call` and scalar indexing patterns in:
   - `/Users/jracine/Development/crs/R/crsiv.R`
   - `/Users/jracine/Development/crs/R/crsivderiv.R`
   - `/Users/jracine/Development/crs/R/crscv.R`
3. Validate with targeted tests + tarball check before next slice.

## Blockers / Decisions Needed

None for planning.  
Implementation is intentionally paused pending your heads-up.
