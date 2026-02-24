# CRS Modernization Plan (R + Non-NOMAD C)

Date: 2026-02-23  
Last Updated: 2026-02-24  
Status: Active execution; major low-risk modernization tranches completed with checkpoint validation.

## Goal

Modernize `crs` to current best-practice R package engineering standards while preserving:

1. Functional behavior.
2. Numerical stability (except accepted optimizer-path drift).
3. Existing user-facing API/contracts.

## Current Status Snapshot (2026-02-24)

1. Checkpoints completed locally: `16+` commits on top of `origin/master` (no push).
2. Completed tranches:
   - `A/R1.1` through `A/R1.10` (low-risk R modernization path).
   - `A/R2.1` (removed `eval(parse(...))` in `stepCV`).
   - `A/R2.2` (removed redundant re-fit work in `add1.lm.cv` path).
   - `B/R1.1` (non-NOMAD C memory hygiene in `gsl_bspline.c`).
3. Validation discipline maintained at each checkpoint:
   - installed-package targeted smokes,
   - tarball-first `R CMD build` and `R CMD check --as-cran`,
   - stable non-regressive check profile (`4 WARNINGs, 2-3 NOTEs`) with no modernization regressions introduced.
4. R-layer forensic status (post A/R1.10):
   - `eval(parse(...))`: `0`
   - string `do.call("...")`: `0`
   - `<<-`: `0`
   - active `1:length(...)`: `0`
   - active `1:ncol(...)`: `0`
   - active `1:NCOL(...)`: `0`
5. Scope guard respected:
   - no edits to NOMAD core source/interface (`src/nomad4_src/**`, `src/snomadr.cpp`, `src/snomadr.h`).

## Accomplished Tasks (Clear Summary)

1. Option-access hardening:
   - migrated fragile `options('x')$x` reads to `getOption(...)` / `isTRUE(...)`.
2. Dynamic evaluation/call cleanup:
   - replaced string `do.call("fn", ...)` dispatch with function references,
   - removed `eval(parse(...))` from `R/stepCV.R`.
3. Scalar safety and control-flow hardening:
   - converted low-risk `1:length`/`1:ncol`/`1:NCOL` patterns to `seq_along`/`seq_len`,
   - converted selected scalar `|`/`&` controls to `||`/`&&`.
4. Super-assignment retirement in R:
   - removed `<<-` state mutation from `console`, `plot.crsiv`, `frscv`, `krscv`, and `np.regression.glp` ridge loop state.
5. Non-NOMAD C hygiene:
   - fixed `gsl_bspline_deriv` allocation cleanup (`gsl_vector_free(quantile_vec)` restored).
6. Project documentation hygiene:
   - each completed tranche logged in the checkpoint section below with concrete `/tmp` artifact paths.

## Next High-ROI Items (Remaining)

1. `A/R2.3`: replace fragile `eval(object$call$data)` handling with robust call/data resolution helpers.
2. `A/R2.4`: reduce repeated call/list rebuilding in `crsiv`/`crsivderiv` iterative paths.
3. `B/R2`: controlled DRY refactors for shared `crsiv`/`crsivderiv` scaffolding after parity checks.

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

No hard blockers at present.  
Execution is active and checkpoint-driven; pending items are listed in `Next High-ROI Items (Remaining)`.

## Checkpoint Log

### 2026-02-23 - A/R1.1 option-access hardening

Scope completed:

1. Replaced direct `options('crs.messages')$crs.messages` reads with scalar-safe `isTRUE(getOption("crs.messages"))` in:
   - `/Users/jracine/Development/crs/R/clsd.R`
   - `/Users/jracine/Development/crs/R/console.R`
   - `/Users/jracine/Development/crs/R/crs.R`
   - `/Users/jracine/Development/crs/R/frscvNOMAD.R`
   - `/Users/jracine/Development/crs/R/krscvNOMAD.R`
   - `/Users/jracine/Development/crs/R/np.regression.glp.R`
2. Modernized startup option reads in:
   - `/Users/jracine/Development/crs/R/zzz.R`
   - `options('np.messages')$np.messages` -> `getOption("np.messages")`
   - `options('np.tree')$np.tree` -> `getOption("np.tree")`
3. Preserved restore semantics where needed:
   - `old.crs.messages <- getOption("crs.messages")` before temporary override.

Validation artifacts:

1. Parse gate:
   - inline run result: `CRS_PARSE_OK`
2. Install and runtime smoke:
   - `/tmp/crs_install_optaccess_20260223.log`
   - `/tmp/crs_smoke_optaccess_20260223.out` (`CRS_SMOKE_OK`)
   - `/tmp/crs_npglpreg_smoke_optaccess_20260223.out` (`NPGLPREG_SMOKE_OK`)
3. Tarball-first:
   - `/tmp/crs_build_optaccess_20260223.log`
   - `/tmp/crs_check_ascran_optaccess_20260223.log` (`Status: 4 WARNINGs, 2 NOTEs`)

Notes:

1. `testthat::test_dir("tests/testthat")` is not authoritative here because tests call unqualified package functions without attaching `crs`; installed-package smoke and tarball checks were used as checkpoint gates.

### 2026-02-23 - A/R1.2 string `do.call` retirement

Scope completed:

1. Replaced string-based call dispatch with function references:
   - `do.call("crs", ...)` -> `do.call(crs, ...)`
   - `do.call("expand.grid", ...)` -> `do.call(base::expand.grid, ...)`
2. Files updated:
   - `/Users/jracine/Development/crs/R/crsiv.R`
   - `/Users/jracine/Development/crs/R/crsivderiv.R`
   - `/Users/jracine/Development/crs/R/np.regression.glp.R`

Validation artifacts:

1. Parse gate:
   - inline run result: `DOCALL_PARSE_OK`
2. Deterministic install/smoke using explicit library:
   - `/tmp/crs_install_docall_20260223.log`
   - `/tmp/crsiv_smoke_docall_20260223.out` (`CRSIV_SMOKE_OK`)
   - `/tmp/crs_npglpreg_smoke_docall_20260223.out` (`NPGLPREG_SMOKE_OK`)
3. Tarball-first:
   - `/tmp/crs_build_docall_20260223.log`
   - `/tmp/crs_check_ascran_docall_20260223.log` (`Status: 4 WARNINGs, 2 NOTEs`)

### 2026-02-23 - A/R1.3 index/loop safety (`seq_len`/`seq_along`)

Scope completed:

1. Replaced selected `1:ncol(...)` / `1:length(...)` loop/index patterns with zero-length-safe forms in low-risk paths:
   - `/Users/jracine/Development/crs/R/crscv.R`
   - `/Users/jracine/Development/crs/R/crs.R`
   - `/Users/jracine/Development/crs/R/crsiv.R`
   - `/Users/jracine/Development/crs/R/crsivderiv.R`
   - `/Users/jracine/Development/crs/R/print.snomadr.R`
   - `/Users/jracine/Development/crs/R/get.option.types.R`
   - `/Users/jracine/Development/crs/R/util.R`
2. Included normalization-index safety:
   - `1:length(norm.stop)` -> `seq_along(norm.stop)` in IV iterative summaries.

Validation artifacts:

1. Parse gate:
   - inline run result: `SEQLEN_PARSE_OK`
2. Deterministic install/smokes (`R CMD INSTALL -l /tmp/crs_lib_seqlen_20260223`):
   - `/tmp/crs_install_seqlen_20260223.log`
   - `/tmp/crs_smoke_seqlen_20260223.out` (`CRS_SMOKE_OK`)
   - `/tmp/crsiv_smoke_seqlen_20260223.out` (`CRSIV_SMOKE_OK`)
   - `/tmp/crsivderiv_smoke_seqlen_20260223.out` (`CRSIVDERIV_SMOKE_OK`)
   - `/tmp/crs_npglpreg_smoke_seqlen_20260223.out` (`NPGLPREG_SMOKE_OK`)
3. Tarball-first:
   - `/tmp/crs_build_seqlen_20260223.log`
   - `/tmp/crs_check_ascran_seqlen_20260223.log` (`Status: 4 WARNINGs, 2 NOTEs`)

### 2026-02-23 - A/R2.1 remove `eval(parse(...))` in `stepCV`

Scope completed:

1. Replaced `eval(parse(...))` formula construction in `R/stepCV.R` with structured formula creation:
   - `stats::reformulate(c(".", scope), env = environment(formula(object)))`
2. Updated both affected callsites in:
   - `/Users/jracine/Development/crs/R/stepCV.R`

Validation artifacts:

1. Parse gate:
   - inline run result: `STEPCV_PARSE_OK`
2. Deterministic install/smoke:
   - `/tmp/crs_install_stepcv_20260223.log`
   - `/tmp/crs_stepcv_smoke_20260223.out` (`STEPCV_SMOKE_OK`, invoked via `crs:::stepCV`)
3. Tarball-first:
   - `/tmp/crs_build_stepcv_20260223.log`
   - `/tmp/crs_check_ascran_stepcv_20260223.log` (`Status: 4 WARNINGs, 2 NOTEs`)

### 2026-02-23 - B/R1.1 non-NOMAD C memory hygiene (`gsl_bspline`)

Scope completed:

1. Fixed leak candidate in derivative wrapper by restoring missing free:
   - `/Users/jracine/Development/crs/src/gsl_bspline.c`
   - `gsl_vector_free(quantile_vec);` in `gsl_bspline_deriv(...)`

Validation artifacts:

1. Deterministic install:
   - `/tmp/crs_install_gslfree_20260223.log`
2. Focused smokes:
   - `/tmp/crs_gslbs_smoke_gslfree_20260223.out` (`GSL_BS_SMOKE_OK`)
   - `/tmp/crs_basic_smoke_gslfree_20260223.out` (`CRS_BASIC_SMOKE_OK`)
3. Tarball-first:
   - `/tmp/crs_build_gslfree_20260223.log`
   - `/tmp/crs_check_ascran_gslfree_20260223.log` (`Status: 4 WARNINGs, 2 NOTEs`)

### 2026-02-23 - A/R1.4 scalar logical operator hardening

Scope completed:

1. Replaced scalar control operators in selected low-risk paths:
   - `&` -> `&&`
   - `|` -> `||`
2. Files updated:
   - `/Users/jracine/Development/crs/R/console.R`
   - `/Users/jracine/Development/crs/R/clsd.R`
   - `/Users/jracine/Development/crs/R/crs.R`

Validation artifacts:

1. Parse gate:
   - inline run result: `SCALAROP_PARSE_OK`
2. Deterministic install/smokes:
   - `/tmp/crs_install_scalarop_20260223.log`
   - `/tmp/crs_clsd_smoke_scalarop_20260223.out` (`CLSD_SMOKE_OK`)
   - `/tmp/crs_basic_smoke_scalarop_20260223.out` (`CRS_BASIC_SMOKE_OK`)
3. Tarball-first:
   - `/tmp/crs_build_scalarop_20260223.log`
   - `/tmp/crs_check_ascran_scalarop_20260223.log` (`Status: 4 WARNINGs, 2 NOTEs`)

### 2026-02-23 - A/R1.5 retire `<<-` in console tab formatting

Scope completed:

1. Removed closure mutation in `toMsg()` tab-fill logic:
   - replaced `tmpPos <<- ...` updates inside `sapply` with explicit local `for` loop accumulation.
2. File updated:
   - `/Users/jracine/Development/crs/R/console.R`

Validation artifacts:

1. Parse gate:
   - inline run result: `CONSOLE_PARSE_OK`
2. Deterministic install/smoke:
   - `/tmp/crs_install_console_20260223.log`
   - `/tmp/crs_console_smoke_20260223.out` (`CONSOLE_SMOKE_OK`)
3. Tarball-first:
   - `/tmp/crs_build_console_20260223.log`
   - `/tmp/crs_check_ascran_console_20260223.log` (`Status: 4 WARNINGs, 2 NOTEs`)

### 2026-02-23 - A/R1.6 retire `<<-` in `plot.crsiv` dot consumption

Scope completed:

1. Removed super-assignment-based state mutation in `plot.crsiv(...)` dot handling:
   - replaced `dots[[name]] <<- NULL` pattern with a local environment (`dot_env`) and explicit consume helpers.
2. Preserved argument-consumption behavior by passing only unconsumed args at each plotting call via `remaining.dots()`.
3. File updated:
   - `/Users/jracine/Development/crs/R/crsiv.R`

Validation artifacts:

1. Deterministic install:
   - `/tmp/crs_install_crsivdots_20260223.log`
2. Focused smoke:
   - `/tmp/crsiv_plot_smoke_crsivdots_20260223.out` (`CRSIV_PLOT_SMOKE_OK`)
   - `/tmp/crsiv_plot_crsivdots_20260223.pdf`
3. Tarball-first:
   - `/tmp/crs_build_crsivdots_20260223.log`
   - `/tmp/crs_check_ascran_crsivdots_20260223.log` (`Status: 4 WARNINGs, 2 NOTEs`)

### 2026-02-23 - A/R1.7 index/loop safety sweep (additional low-risk paths)

Scope completed:

1. Replaced additional `1:length(...)`, `1:ncol(...)`, and `1:NCOL(...)` loop/index patterns with `seq_along(...)` / `seq_len(...)` in low-risk paths.
2. Files updated:
   - `/Users/jracine/Development/crs/R/np.regression.glp.R`
   - `/Users/jracine/Development/crs/R/clsd.R`
   - `/Users/jracine/Development/crs/R/snomadr.R`
   - `/Users/jracine/Development/crs/R/crssigtest.R`

Validation artifacts:

1. Deterministic install:
   - `/tmp/crs_install_indexsafe2_20260223.log`
2. Focused smoke (touched surfaces):
   - `/tmp/crs_indexsafe2_smoke_20260223.out` (`INDEXSAFE2_SMOKE_OK`)
3. Tarball-first:
   - `/tmp/crs_build_indexsafe2_20260223.log`
   - `/tmp/crs_check_ascran_indexsafe2_20260223.log` (`Status: 4 WARNINGs, 2 NOTEs`)

### 2026-02-23 - A/R1.8 retire `<<-` in `frscv`/`krscv` console updates

Scope completed:

1. Removed `console <<- ...` closure mutation in cross-validation objective callbacks for:
   - `/Users/jracine/Development/crs/R/frscv.R`
   - `/Users/jracine/Development/crs/R/krscv.R`
2. Replaced with explicit local state container (`console.state$console`) that is mutated directly without super-assignment.
3. Kept user-facing progress behavior unchanged.

Validation artifacts:

1. Deterministic install:
   - `/tmp/crs_install_consoleenv_20260223.log`
2. Focused smoke:
   - `/tmp/crs_frs_krs_smoke_consoleenv_20260223.out` (`FRS_KRS_CONSOLEENV_SMOKE_OK`)
3. Tarball-first:
   - `/tmp/crs_build_consoleenv_20260223.log`
   - `/tmp/crs_check_ascran_consoleenv_20260223.log` (`Status: 4 WARNINGs, 2 NOTEs`)

### 2026-02-23 - A/R1.9 retire last R-layer `<<-` usage (`np.regression.glp`)

Scope completed:

1. Replaced ridge-loop super-assignment state in `glpregEst(...)` with explicit local environment state:
   - `ridge`, `ridge.lc`, and `doridge` now managed via `ridge.state$...`.
2. Updated related loops to use `seq_len(...)` where appropriate in the same block.
3. File updated:
   - `/Users/jracine/Development/crs/R/np.regression.glp.R`
4. Post-change audit:
   - `rg -n "<<-" /Users/jracine/Development/crs/R` returns no matches.

Validation artifacts:

1. Deterministic install:
   - `/tmp/crs_install_ridgestate_20260223.log`
2. Focused smoke:
   - `/tmp/crs_npglp_smoke_ridgestate_20260223.out` (`RIDGESTATE_SMOKE_OK`)
3. Tarball-first:
   - `/tmp/crs_build_ridgestate_20260223.log`
   - `/tmp/crs_check_ascran_ridgestate_20260223.log` (`Status: 4 WARNINGs, 2 NOTEs`)

### 2026-02-23 - A/R1.10 residual index-range sweep (`crs`/`spline`/`util`)

Scope completed:

1. Replaced remaining active `1:NCOL(...)` / `1:length(...)` patterns with `seq_len(...)` / `seq_along(...)` in:
   - `/Users/jracine/Development/crs/R/crs.R`
   - `/Users/jracine/Development/crs/R/spline.R`
   - `/Users/jracine/Development/crs/R/util.R`
2. Post-change audit:
   - active `1:NCOL(` / `1:length(` matches removed from code paths (remaining matches are comments in `clsd.R` only).

Validation artifacts:

1. Deterministic install:
   - `/tmp/crs_install_lastindex_20260223.log`
2. Focused smoke:
   - `/tmp/crs_lastindex_smoke_20260223.out` (`LASTINDEX_SMOKE_OK`)
3. Tarball-first:
   - `/tmp/crs_build_lastindex_20260223.log`
   - `/tmp/crs_check_ascran_lastindex_20260223.log` (`Status: 4 WARNINGs, 2 NOTEs`)

### 2026-02-23 - Forensic Snapshot (post A/R1.10)

Current R-layer pattern counts:

1. `eval(parse(...))`: `0`
2. `do.call("<string>", ...)`: `0`
3. `<<-`: `0`
4. `1:length(...)`: `0` in active code (`2` comment-only occurrences in `clsd.R`)
5. `1:ncol(...)`: `0`
6. `1:NCOL(...)`: `0` in active code
7. `.C(` callsites in `R/`: `0`
8. `.Call(` callsites in `R/`: `0`

### 2026-02-24 - A/R2.2 remove redundant `lm()` re-fit in `add1.lm.cv`

Scope completed:

1. Removed duplicate per-candidate `lm(y~X)` fit in `add1.lm.cv(...)` and replaced with a single helper path that computes LOO CV from matrix inputs.
2. Added internal helper `.loocv_from_xy(X, y)` in:
   - `/Users/jracine/Development/crs/R/stepCV.R`
3. Implementation details:
   - uses `.lm.fit` on `cbind("(Intercept)" = 1, X)` (no formula re-parse overhead),
   - computes hat diagonals via `.Call("crs_hat_diag", ...)` when available,
   - includes QR fallback with rank-aware `Q[, seq_len(rank), drop = FALSE]` to preserve numerical parity in rank-deficient cases.

Validation artifacts:

1. Objective/parity probe (old vs new):
   - `/tmp/crs_stepcv_parity_probe_20260224.R`
   - `/tmp/crs_stepcv_parity_probe_20260224.out`
   - result: `MAX_ABS_CV_DIFF=0`, `ALL_EQUAL_CV=TRUE`
2. Focused micro-benchmark (`add1.lm.cv` loop):
   - `/tmp/crs_stepcv_add1_bench_20260224.R`
   - `/tmp/crs_stepcv_add1_bench_20260224.out`
   - result: `OLD_ELAPSED=0.310000`, `NEW_ELAPSED=0.219000`, `DELTA_PCT=-29.35`
3. Deterministic install/smoke:
   - `/tmp/crs_install_stepcvperf_20260224.log`
   - `/tmp/crs_stepcvperf_smoke_20260224.out` (`STEPCVPERF_SMOKE_OK`)
4. Tarball-first:
   - `/tmp/crs_build_stepcvperf_20260224.log`
   - `/tmp/crs_check_ascran_stepcvperf_20260224.log` (`Status: 4 WARNINGs, 3 NOTEs`)
