# CRS Modernization Plan (R + Non-NOMAD C)

Date: 2026-02-23  
Last Updated: 2026-02-24  
Status: Active execution; low-risk modernization scope is effectively complete with checkpoint validation, and remaining items are selective/high-risk refactors.

## Goal

Modernize `crs` to current best-practice R package engineering standards while preserving:

1. Functional behavior.
2. Numerical stability (except accepted optimizer-path drift).
3. Existing user-facing API/contracts.

## Current Status Snapshot (2026-02-24)

1. Checkpoints completed locally: `75+` commits on top of `origin/master` (no push).
2. Completed tranches:
   - `A/R1.1` through `A/R1.10` (low-risk R modernization path).
   - `A/R2.1` (removed `eval(parse(...))` in `stepCV`).
   - `A/R2.2` (removed redundant re-fit work in `add1.lm.cv` path).
   - `A/R2.3` (hardened `crs.sigtest` call-data resolution).
   - `A/R2.4` (reduced repeated call/list rebuilding in `crsiv`/`crsivderiv` loops).
   - `A/R2.5` (removed remaining direct `eval(...)` callsites in `R/` via parity-safe `.crs_eval_call` helper and new regression tests).
   - `A/R1.11` (scalar stopping-rule branch hardening in IV iteration loops).
   - `A/R1.12` (additional `seq_len`/`seq.int` safety sweep in spline utilities).
   - `A/R1.13` (broader loop-header safety sweep across CV/sigtest and helper paths).
   - `A/R1.14` (full `npglpreg` loop-header safety sweep).
   - `A/R1.15` through `A/R1.52` (scalar control-flow cleanup, script hygiene, clamp/vectorization follow-ups, remaining legacy `do.call` cleanup in matrix-construction paths, `vapply` numeric-column sweeps, boolean control simplification, centralized RNG seed state handling, scalar-loop index correctness fixes in `npglpreg` warning/index paths, registered-symbol `.Call` hygiene, robust option-state guards for `crsiv` / `crsivderiv` including unset-option handling, eval-helper namespace-resolution hardening, recursive source-artifact cleanup/build hygiene hardening, stepCV loop-index safety hardening for empty-scope paths, snomadr argument-validation hardening, stepCV call-evaluation/helper cleanup, snomadr single-argument `...` guard hardening, formal-matching parity for defaults/ellipsis, non-NOMAD C forensic comment cleanup, residual scalar bitwise-operator cleanup in control paths, `plot.clsd` scalar logical follow-up alignment, invariant-hoisting in `W.glp` polynomial index construction, explicit `TRUE`/`FALSE` constant normalization in active call arguments, duplicate helper-definition consolidation for `scale_robust`/`is.fullrank`, `on.exit(..., add=TRUE)` stacking safety for `npglpreg` message-state restoration, sequence-construction hardening in GLP matrix expansion, additional zero-dimension-safe sequence/rep-index hardening in `glp.model.matrix`, residual `1:n` index-slice hardening in GLP/NPGLP bandwidth internals, `frscv/krscv` family index-slice normalization to `seq_len` in R wrappers, final remaining executable `1:num.*`/`[1:n]` slice hardening in `frscv`/`krscv`/`spline`, closure of the last executable `1:num.*` slice sites across these wrappers with static-audit confirmation, parity-safe `seq.int(2L, ...)` normalization for remaining guarded loop headers in active non-NOMAD R paths, and guarded `iterate.max` loop-header hardening in `crsiv`/`crsivderiv` with deterministic convergence defaults).
   - `B/R2` (shared IV scaffolding helpers for dots/call assembly in `crsiv` and `crsivderiv`).
   - `B/R1.1` (non-NOMAD C memory hygiene in `gsl_bspline.c`).
   - `B/R1.2` and `B/R1.3` (native-interface and matrix-kernel regression test expansion).
   - test harness stabilization (`tests/testthat/setup-load-crs.R`) and deprecation cleanup (`test-crsivderiv.R`).
   - `C/R3.1` and `C/R3.2` partial migrations: `uniquecombs()` and `gsl.bs` native paths moved from `.C` to `.Call`.
   - `C/R3.3` parity tranche: explicit `.Call` vs legacy `.C` regression tests for `bs.des`.
3. Validation discipline maintained at each checkpoint:
   - installed-package targeted smokes,
   - tarball-first `R CMD build` and `R CMD check --as-cran`,
   - stable non-regressive check profile (`4-5 WARNINGs, 0-4 NOTEs`) with no modernization regressions introduced,
   - latest post-tranche gate in this session: `Status: 5 WARNINGs` (the intermittent `unable to verify current time` NOTE did not appear in that run).
4. R-layer forensic status (post A/R2.5 completion):
   - `eval(parse(...))`: `0`
   - `eval(...)`: `0`
   - string `do.call("...")`: `0`
   - `<<-`: `0`
   - active `1:length(...)`: `0`
   - active `1:ncol(...)`: `0`
   - active `1:NCOL(...)`: `0`
   - scalar `ifelse(is.null(...))` / `ifelse(is.finite(...))` / `return(ifelse(...))`: `0`
   - total `ifelse(...)` uses in `R/`: `0`
   - `.C(` callsites in `R/`: `0`
   - `.Call(` callsites in `R/`: `10`
   - active `for (... in 2:...)` loop headers in `R/`: `0`
5. Scope guard respected:
   - no edits to NOMAD core source/interface (`src/nomad4_src/**`, `src/snomadr.cpp`, `src/snomadr.h`).

## Session Addendum (2026-02-24, latest)

1. Completed and committed tranche checkpoints:
   - `e69fbee` `modernize(r): normalize guarded 2:n loop headers with seq.int`
   - `66b91af` `modernize(iv): guard iterate.max loops and normalize loop headers`
2. Static forensic sweeps after these checkpoints:
   - no executable `for (... in 2:...)` loop headers remain in `R/`,
   - legacy index-pattern scan remains clean; remaining `1:length(...)` hits are comment-only in `/Users/jracine/Development/crs/R/clsd.R`.
3. Validation artifacts for the latest tranche gates:
   - `/tmp/crs_parse_loop_seqint_guarded_20260224.out`
   - `/tmp/crs_install_loop_seqint_guarded_20260224.log`
   - `/tmp/crs_test_loop_seqint_guarded_targeted_20260224.out`
   - `/tmp/crs_test_loop_seqint_guarded_full_20260224.out`
   - `/tmp/crs_build_loop_seqint_guarded_20260224.log`
   - `/tmp/crs_check_loop_seqint_guarded_20260224.log`
   - `/tmp/crs_parse_iterate_seq_guard_clean_20260224.out`
   - `/tmp/crs_install_iterate_seq_guard_20260224.log`
   - `/tmp/crs_test_iterate_seq_guard_targeted_20260224.out`
   - `/tmp/crs_test_iterate_seq_guard_full_20260224.out`
   - `/tmp/crs_build_iterate_seq_guard_20260224.log`
   - `/tmp/crs_check_iterate_seq_guard_20260224.log`

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

1. `C/R3` (optional/deferred): deeper native-interface API cleanup and algorithmic rewrites only after extended parity/performance signoff.

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

## Baseline Audit Snapshot (`crs`, 2026-02-23 initial)

Static inventory (initial baseline):

1. `.C(` callsites in `R/`: `3`
2. `.Call(` callsites in `R/`: `6`
3. `eval(parse(...))` in `R/`: `2` (both in `R/stepCV.R`)
4. `eval(` in `R/`: `7`
5. `do.call("<string>", ...)` in `R/`: `13`
6. `<<-` in `R/`: `15`
7. `ifelse(` in `R/`: `58`
8. `1:length(...)` in `R/`: `50`

Current inventory (2026-02-24):

1. `.C(` callsites in `R/`: `0`
2. `.Call(` callsites in `R/`: `10`
3. `eval(parse(...))` in `R/`: `0`
4. `eval(` in `R/`: `0`
5. string `do.call("<name>", ...)` in `R/`: `0`
6. `<<-` in `R/`: `0`
7. `ifelse(` in `R/`: `0`

Resolved hotspot history (non-NOMAD-heavy):

1. `/Users/jracine/Development/crs/R/spline.R` (`ifelse` density) resolved.
2. `/Users/jracine/Development/crs/R/crs.R` (`eval(object$call$data)`, indexing patterns) resolved for direct `eval(...)`.
3. `/Users/jracine/Development/crs/R/crsiv.R` and `/Users/jracine/Development/crs/R/crsivderiv.R` (string `do.call`, repeated loop scaffolding) resolved.
4. `/Users/jracine/Development/crs/R/stepCV.R` (`eval(parse(...))`, duplicate model fitting in CV loop) resolved.
5. `/Users/jracine/Development/crs/R/console.R`, `/Users/jracine/Development/crs/R/frscv.R`, `/Users/jracine/Development/crs/R/krscv.R`, `/Users/jracine/Development/crs/R/np.regression.glp.R` (`<<-` usage) resolved.

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
7. `.C(` callsites in `R/`: `3`
8. `.Call(` callsites in `R/`: `7`

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

### 2026-02-24 - A/R2.3 harden `crs.sigtest` call-data resolution

Scope completed:

1. Replaced direct `eval(object$call$data)` usage in `crs.sigtest(...)` with a local resolver helper:
   - `resolve.call.data(object)` now resolves `object$call$data` across robust environment candidates:
     - `attr(object$terms, ".Environment")`
     - `environment(object$formula)`
     - `parent.frame()`
2. Resolved data object once per invocation (`model.data`) and reused it in the variable loop.
3. Updated callsites in:
   - `/Users/jracine/Development/crs/R/crs.R`

Validation artifacts:

1. Parse gate:
   - inline run result: `CRS_PARSE_OK`
2. Deterministic install:
   - `/tmp/crs_install_sigtest_20260224.log`
3. Focused smoke:
   - `/tmp/crs_sigtest_smoke_20260224.out` (`CRS_SIGTEST_SMOKE_OK`)
4. Regression smoke for detached data-symbol scope (local function fit, later `crs.sigtest`):
   - `/tmp/crs_sigtest_envresolve_20260224.R`
   - `/tmp/crs_sigtest_envresolve_20260224.out` (`CRS_SIGTEST_ENV_RESOLVE_OK`)
5. Tarball-first:
   - `/tmp/crs_build_sigtestdata_20260224.log`
   - `/tmp/crs_check_ascran_sigtestdata_20260224.log` (`Status: 4 WARNINGs, 2 NOTEs`)

### 2026-02-24 - A/R2.4 reduce repeated call/list rebuilding in `crsiv*`

Scope completed:

1. Added local `fit.crs(...)` helpers in:
   - `/Users/jracine/Development/crs/R/crsiv.R`
   - `/Users/jracine/Development/crs/R/crsivderiv.R`
2. Replaced repeated `do.call(crs, c(list(...), dots.*))` blocks in pre-loop and iterative paths with helper calls that centralize shared arguments (`opts`, `data`, `display.*`, optional warm-start args).
3. Removed redundant mid-function re-creation of `dots.loop` and reused pre-sanitized loop dots to avoid repeated rebuild overhead and duplicate-arg risk.

Validation artifacts:

1. Parse gate:
   - inline run result: `CRSIV_PARSE_OK`
2. Deterministic install:
   - `/tmp/crs_install_crsivref_20260224.log`
3. Focused runtime smoke (`crsiv` + `crsivderiv` with `cv="none"`, `kernel=FALSE`, small `iterate.max`):
   - `/tmp/crsiv_refactor_smoke_20260224.R`
   - `/tmp/crsiv_refactor_smoke_20260224.out` (`CRSIV_REFACTOR_SMOKE_OK`)
4. Tarball-first:
   - `/tmp/crs_build_crsivref_20260224.log`
   - `/tmp/crs_check_ascran_crsivref_20260224.log` (`Status: 4 WARNINGs, 3 NOTEs`)

### 2026-02-24 - B/R2 shared `crsiv*` scaffolding helpers

Scope completed:

1. Added shared internal helper functions in:
   - `/Users/jracine/Development/crs/R/util.R`
   - `.crsiv_prepare_dot_args(...)` for consistent `...` sanitization across pre-loop and loop calls.
   - `.crsiv_fit_crs(...)` for standardized `crs` call assembly with optional warm-start args.
2. Updated both IV entry points to use shared helpers:
   - `/Users/jracine/Development/crs/R/crsiv.R`
   - `/Users/jracine/Development/crs/R/crsivderiv.R`
3. Kept behavior contracts stable while removing duplicated local scaffolding blocks.

Validation artifacts:

1. Parse gate:
   - inline run result: `CRSIV_B2_PARSE_OK`
2. Deterministic install:
   - `/tmp/crs_install_crsivb2_20260224.log`
3. Focused runtime smoke (`crsiv` + `crsivderiv`):
   - `/tmp/crsiv_b2_smoke_20260224.out` (`CRSIV_REFACTOR_SMOKE_OK`)
4. Tarball-first:
   - `/tmp/crs_build_crsivb2_20260224.log`
   - `/tmp/crs_check_ascran_crsivb2_20260224.log` (`Status: 4 WARNINGs, 2 NOTEs`)

### 2026-02-24 - A/R1.11 scalar stopping-rule branch hardening (`crsiv*`)

Scope completed:

1. Replaced scalar `ifelse(...)` stopping-rule updates in iterative IV loops with scalar `if (...) ... else ...` branches and single-pass raw score computation:
   - `/Users/jracine/Development/crs/R/crsiv.R`
   - `/Users/jracine/Development/crs/R/crsivderiv.R`
2. Reduced duplicate per-iteration arithmetic by computing `norm.raw` once before branch application.

Validation artifacts:

1. Parse gate:
   - inline run result: `CRSIV_NORM_PARSE_OK`
2. Deterministic install:
   - `/tmp/crs_install_normbranch_20260224.log`
3. Focused runtime smoke (`crsiv` + `crsivderiv`):
   - `/tmp/crsiv_normbranch_smoke_20260224.R`
   - `/tmp/crsiv_normbranch_smoke_20260224.out` (`CRSIV_NORMBRANCH_SMOKE_OK`)
4. Tarball-first:
   - `/tmp/crs_build_normbranch_20260224.log`
   - `/tmp/crs_check_ascran_normbranch_20260224.log` (`Status: 4 WARNINGs, 3 NOTEs`)

### 2026-02-24 - A/R1.12 index-range safety sweep (`spline`/`util`)

Scope completed:

1. Replaced remaining targeted `for (... in 1:...)` loop headers with zero-length-safe forms in:
   - `/Users/jracine/Development/crs/R/spline.R`
   - `/Users/jracine/Development/crs/R/util.R`
2. Updated guarded upper-range loop form:
   - `for(i in 2:num.z)` -> `for(i in seq.int(2, num.z))` under `if(num.z > 1)`.

Validation artifacts:

1. Parse gate:
   - inline run result: `SPLINE_SEQ_PARSE_OK`
2. Deterministic install:
   - `/tmp/crs_install_splineseq_20260224.log`
3. Focused spline-path runtime smoke (kernel/non-kernel + eval data + derivative path):
   - `/tmp/crs_spline_seq_smoke_20260224.R`
   - `/tmp/crs_spline_seq_smoke_20260224.out` (`CRS_SPLINE_SEQ_SMOKE_OK`)
4. Tarball-first:
   - `/tmp/crs_build_splineseq_20260224.log`
   - `/tmp/crs_check_ascran_splineseq_20260224.log` (`Status: 4 WARNINGs, 3 NOTEs`)

### 2026-02-24 - A/R1.13 broader loop-header safety sweep

Scope completed:

1. Replaced additional active `for (... in 1:...)` loop headers with zero-length-safe forms in:
   - `/Users/jracine/Development/crs/R/crs.R`
   - `/Users/jracine/Development/crs/R/clsd.R`
   - `/Users/jracine/Development/crs/R/snomadr.R`
   - `/Users/jracine/Development/crs/R/frscv.R`
   - `/Users/jracine/Development/crs/R/glp.model.matrix.R`
   - `/Users/jracine/Development/crs/R/krscv.R`
   - `/Users/jracine/Development/crs/R/krscvNOMAD.R`
   - `/Users/jracine/Development/crs/R/crssigtest.R`
2. Kept scope strictly R-layer and non-invasive (loop-header normalization only; no NOMAD core edits).

Validation artifacts:

1. Parse gate:
   - `/tmp/crs_parse_loopsweep_20260224.log` (`LOOP_SWEEP_PARSE_OK`)
2. Deterministic install:
   - `/tmp/crs_install_loopsweep2_20260224.log`
3. Focused runtime smoke (exhaustive/non-kernel, exhaustive/kernel, NOMAD-kernel, GLP, and `crssigtest` reorder path):
   - `/tmp/crs_loop_sweep_smoke_20260224.R`
   - `/tmp/crs_loop_sweep_smoke_20260224.out` (`CRS_LOOP_SWEEP_SMOKE_OK`)
4. Tarball-first:
   - `/tmp/crs_build_loopsweep_20260224.log`
   - `/tmp/crs_check_ascran_loopsweep_20260224.log` (`Status: 4 WARNINGs, 3 NOTEs`)

### 2026-02-24 - A/R1.14 `npglpreg` loop-header safety sweep

Scope completed:

1. Converted remaining active `for (... in 1:...)` loop headers in:
   - `/Users/jracine/Development/crs/R/np.regression.glp.R`
   to zero-length-safe `seq_len(...)` forms.
2. Post-change R-layer scan now shows no active `for(... in 1:...)` loop headers (remaining match is comment-only in `glp.model.matrix.R`).

Validation artifacts:

1. Parse gate:
   - `/tmp/crs_parse_npglp_loopsweep_20260224.log` (`NPGLP_LOOP_PARSE_OK`)
2. Deterministic install:
   - `/tmp/crs_install_npglp_loopsweep_20260224.log`
3. Focused `npglpreg` runtime smoke (`cv.ls` and `cv.aic`):
   - `/tmp/crs_npglp_loopsweep_smoke_20260224.R`
   - `/tmp/crs_npglp_loopsweep_smoke_20260224.out` (`NPGLP_LOOP_SWEEP_SMOKE_OK`)
4. Tarball-first:
   - `/tmp/crs_build_npglp_loopsweep_20260224.log`
   - `/tmp/crs_check_ascran_npglp_loopsweep_20260224.log` (`Status: 4 WARNINGs, 2 NOTEs`)

### 2026-02-24 - A/R1.15 indexer safety + scalar logical hardening

Scope completed:

1. Replaced remaining active `sapply(1:...)` indexers with zero-length-safe `sapply(seq_len(...), ...)` in:
   - `/Users/jracine/Development/crs/R/np.regression.glp.R`
   - `/Users/jracine/Development/crs/R/spline.R`
2. Hardened scalar branch logic to avoid vectorized-or semantics in:
   - `/Users/jracine/Development/crs/R/krscvNOMAD.R`
   - `complexity=="degree"|complexity=="knots"` -> `complexity=="degree" || complexity=="knots"`.
3. Scope remains non-invasive: no NOMAD core code changes.

Validation artifacts:

1. Parse gate:
   - `/tmp/crs_parse_seqapply_scalarlogic_20260224.log` (`SEQAPPLY_SCALARLOGIC_PARSE_OK`)
2. Deterministic install:
   - `/tmp/crs_install_seqapply_scalarlogic2_20260224.log`
3. Focused runtime smoke (GLP + spline/derivative + `krscvNOMAD` branch path):
   - `/tmp/crs_seqapply_scalarlogic_smoke_20260224.R`
   - `/tmp/crs_seqapply_scalarlogic_smoke2_20260224.out` (`CRS_SEQAPPLY_SCALARLOGIC_SMOKE_OK`)
4. Tarball-first:
   - `/tmp/crs_build_seqapply_scalarlogic_20260224.log`
   - `/tmp/crs_check_ascran_seqapply_scalarlogic_20260224.log` (`Status: 4 WARNINGs, 2 NOTEs`)

### 2026-02-24 - A/R1.16 scalar-control `ifelse` cleanup sweep

Scope completed:

1. Replaced clearly scalar `ifelse(...)` control/dispatch patterns with scalar `if (...) ... else ...` branches in:
   - `/Users/jracine/Development/crs/R/np.regression.glp.R`
   - `/Users/jracine/Development/crs/R/spline.R`
   - `/Users/jracine/Development/crs/R/crs.R`
2. Replaced scalar logical `&` with `&&` in bound-check guard and simplified scalar integer flag assignment in:
   - `/Users/jracine/Development/crs/R/gsl_bspline.R`
3. Reduced repeated scalar-expression evaluation for NOMAD restart count forwarding by introducing one local scalar:
   - `nmulti.nomad <- if(nmulti == 1) 0 else nmulti` in `glpcvNOMAD`.
4. Scope remains non-invasive: no NOMAD core source edits.

Validation artifacts:

1. Parse gate:
   - `/tmp/crs_parse_scalar_ifelse_20260224.log` (`SCALAR_IFELSE_PARSE_OK`)
2. Deterministic install:
   - `/tmp/crs_install_scalar_ifelse_20260224.log`
3. Focused runtime smoke:
   - `/tmp/crs_npglp_scalar_ifelse_smoke_20260224.out` (`NPGLP_SMOKE_OK`)
   - `/tmp/crs_spline_scalar_ifelse_smoke_20260224.out` (`SPLINE_SMOKE_OK`)
   - `/tmp/crs_gsl_scalar_ifelse_smoke_20260224.out` (`GSL_SMOKE_OK`)
4. Tarball-first:
   - `/tmp/crs_build_scalar_ifelse_20260224.log`
   - `/tmp/crs_check_ascran_scalar_ifelse_20260224.log` (`Status: 4 WARNINGs, 3 NOTEs`)

### 2026-02-24 - A/R1.17 scalar-label/control cleanup follow-up

Scope completed:

1. Removed remaining scalar `ifelse(is.finite(...), ...)` clamp in:
   - `/Users/jracine/Development/crs/R/np.regression.glp.R`
   using scalar branch `if(!is.finite(fv)) fv <- cv.maxPenalty`.
2. Replaced scalar plot label/title `ifelse(...)` calls in:
   - `/Users/jracine/Development/crs/R/crs.R`
   with scalar `if (...) ... else ...` expressions.
3. Corrected a temporary parse issue in the 3D plot title expression and revalidated package build/check.

Validation artifacts:

1. Parse gate:
   - `/tmp/crs_parse_scalar_ifelse_followup_20260224.log` (`SCALAR_IFELSE_FOLLOWUP_PARSE_OK`)
2. Deterministic install:
   - `/tmp/crs_install_scalar_ifelse_followup_20260224.log`
3. Focused runtime smoke:
   - `/tmp/crs_npglp_scalar_ifelse_followup_smoke_20260224.out` (`NPGLP_SMOKE_OK`)
   - `/tmp/crs_loop_scalar_ifelse_followup_smoke_20260224.out` (`LOOP_SMOKE_OK`)
4. Tarball-first:
   - `/tmp/crs_build_scalar_ifelse_followup_20260224.log`
   - `/tmp/crs_check_ascran_scalar_ifelse_followup_20260224.log` (`Status: 4 WARNINGs, 3 NOTEs`)

### 2026-02-24 - A/R1.18 residual scalar `ifelse` cleanup (`spline`/`npglpreg`)

Scope completed:

1. Replaced remaining scalar `ifelse` constructs in:
   - `/Users/jracine/Development/crs/R/spline.R`
   - `/Users/jracine/Development/crs/R/np.regression.glp.R`
2. Specific updates:
   - scalar confidence-band `se` column assembly in `spline` now uses scalar `if (...)` control flow.
   - final additive derivative index start in `spline` now uses scalar `if (...)` branch.
   - gradient plot labels in `npglpreg` now use scalar `if (...) ... else ...` expressions.
3. Scope remains behavior-preserving and wrapper-only (no NOMAD core edits).

Validation artifacts:

1. Parse gate:
   - `/tmp/crs_parse_scalar_ifelse_finalsweep_20260224.log` (`SCALAR_IFELSE_FINAL_SWEEP_PARSE_OK`)
2. Deterministic install:
   - `/tmp/crs_install_scalar_ifelse_finalsweep_20260224.log`
3. Focused runtime smoke:
   - `/tmp/crs_npglp_scalar_ifelse_finalsweep_smoke_20260224.out` (`NPGLP_SMOKE_OK`)
   - `/tmp/crs_spline_scalar_ifelse_finalsweep_smoke_20260224.out` (`SPLINE_SMOKE_OK`)
4. Tarball-first:
   - `/tmp/crs_build_scalar_ifelse_finalsweep_20260224.log`
   - `/tmp/crs_check_ascran_scalar_ifelse_finalsweep_20260224.log` (`Status: 4 WARNINGs, 3 NOTEs`)

### 2026-02-24 - A/R1.19 script hygiene (`runcrs`/`dontruncrs`)

Scope completed:

1. Hardened backup cleanup command in:
   - `/Users/jracine/Development/crs/man/runcrs`
   - `/Users/jracine/Development/crs/man/dontruncrs`
2. Changed:
   - `rm *.Rdbak` -> `rm -f -- ./*.Rdbak`
   to avoid option-interpretation edge cases and satisfy `shellcheck` guidance.

Validation artifacts:

1. Script lint/syntax:
   - `shellcheck /Users/jracine/Development/crs/man/runcrs /Users/jracine/Development/crs/man/dontruncrs` (clean)
   - `bash -n /Users/jracine/Development/crs/man/runcrs` and `bash -n /Users/jracine/Development/crs/man/dontruncrs` (`BASH_N_OK`)
2. Tarball-first:
   - `/tmp/crs_build_script_hygiene_20260224.log`
   - `/tmp/crs_check_ascran_script_hygiene_20260224.log` (`Status: 4 WARNINGs, 3 NOTEs`)

### 2026-02-24 - A/R1.20 clamp/vectorization normalization sweep

Scope completed:

1. Replaced equivalent clamp/select `ifelse(...)` patterns with `pmin(...)`/`pmax(...)` or direct logical expressions in:
   - `/Users/jracine/Development/crs/R/np.regression.glp.R`
   - `/Users/jracine/Development/crs/R/krscv.R`
   - `/Users/jracine/Development/crs/R/frscvNOMAD.R`
   - `/Users/jracine/Development/crs/R/krscvNOMAD.R`
   - `/Users/jracine/Development/crs/R/crssigtest.R`
   - `/Users/jracine/Development/crs/R/spline.R`
2. Normalized degree/lambda upper/lower bound clamps and hat-value caps via vectorized intrinsics (`pmin`/`pmax`) for clearer intent and lower overhead.
3. Replaced boolean-selection `ifelse` with direct logical algebra where appropriate (`deriv.ind.vec <- deriv.ind.vec & prune.index`).
4. R-layer `ifelse(...)` inventory reduced from `29` to `2` (remaining intentional vectorized cases in `util.R` and `spline.R` helper math).

Validation artifacts:

1. Parse gate:
   - `/tmp/crs_parse_clamp_sweep_20260224.log` (`CLAMP_SWEEP_PARSE_OK`)
2. Deterministic install:
   - `/tmp/crs_install_clamp_sweep_20260224.log`
3. Focused runtime smoke:
   - `/tmp/crs_npglp_clamp_sweep_smoke_20260224.out` (`NPGLP_SMOKE_OK`)
   - `/tmp/crs_spline_clamp_sweep_smoke_20260224.out` (`SPLINE_SMOKE_OK`)
   - `/tmp/crs_loop_clamp_sweep_smoke_20260224.out` (`LOOP_SMOKE_OK`)
   - `/tmp/crs_clamp_sweep_nomad_smoke_20260224.out` (`NOMAD_CV_SMOKE_OK`)
4. Tarball-first:
   - `/tmp/crs_build_clamp_sweep_20260224.log`
   - `/tmp/crs_check_ascran_clamp_sweep_20260224.log` (`Status: 4 WARNINGs, 3 NOTEs`)

### 2026-02-24 - Periodic heavy-example gate (`runcrs`)

Scope completed:

1. Executed full example stress gate using package-maintained toggles:
   - `/Users/jracine/Development/crs/man/runcrs`
   - `R CMD build` + `R CMD check --as-cran` from `/Users/jracine/Development`
   - `/Users/jracine/Development/crs/man/dontruncrs`
2. Verified `dontruncrs` restoration completed successfully after check run.

Validation artifacts:

1. Toggle/build/check logs:
   - `/tmp/crs_runcrs_toggle_on_20260224.log` (`RUNCRS_ON_OK`)
   - `/tmp/crs_build_runcrs_gate_20260224.log` (`BUILD_OK`)
   - `/tmp/crs_check_ascran_runcrs_gate_20260224.log` (`Status: 4 WARNINGs, 4 NOTEs`)
   - `/tmp/crs_runcrs_toggle_off_20260224.log` (`RUNCRS_OFF_OK`)
2. Example outputs:
   - `/Users/jracine/Development/crs.Rcheck/crs-Ex.R`
   - `/Users/jracine/Development/crs.Rcheck/crs-Ex.Rout`
   - `/Users/jracine/Development/crs.Rcheck/00check.log`
3. Note profile change vs standard gate:
   - additional NOTE from long examples timing (`* checking examples ... [11m/11m] NOTE`) under full-example mode.

### 2026-02-24 - B/R1.2 native-interface regression tests

Scope completed:

1. Added focused non-NOMAD native-interface regression tests in:
   - `/Users/jracine/Development/crs/tests/testthat/test-native-interfaces.R`
2. Coverage added for:
   - `RuniqueCombs` path via `uniquecombs()` row-index reconstruction invariants.
   - `gsl_bspline`/`gsl_bspline_deriv` paths via `gsl.bs()` fit/predict dimension invariants.
   - `crs_hat_diag` path via `crs:::hat.from.lm.fit()` parity against QR-based `hat(...)`.
3. Test file uses explicit namespace qualification (`crs::` / `crs:::`) to support direct `test_dir(...)` execution.

Validation artifacts:

1. Parse + local test execution:
   - `/tmp/crs_parse_native_tests_20260224.log` (`NATIVE_TEST_PARSE_OK`)
   - `/tmp/crs_test_native_interfaces_20260224.out` (`NATIVE_TESTS_OK`)
2. Deterministic install and smoke:
   - `/tmp/crs_install_native_interfaces_20260224.log`
   - `/tmp/crs_loop_native_interfaces_smoke_20260224.out` (`LOOP_SMOKE_OK`)
3. Tarball-first:
   - `/tmp/crs_build_native_interfaces_20260224.log`
   - `/tmp/crs_check_ascran_native_interfaces_20260224.log` (`Status: 4 WARNINGs, 3 NOTEs`)

### 2026-02-24 - Test harness stabilization for direct `test_dir` runs

Scope completed:

1. Added test setup bootstrap:
   - `/Users/jracine/Development/crs/tests/testthat/setup-load-crs.R`
2. Purpose:
   - ensure direct local `testthat::test_dir('tests/testthat')` runs have package symbols loaded, matching expected developer workflow outside `R CMD check`.
3. Result:
   - full local test suite now runs successfully in direct mode.

Validation artifacts:

1. Full-suite direct test run:
   - `/tmp/crs_test_fullsuite_20260224.out` (`FULL_TESTS_OK`)
2. Tarball-first:
   - `/tmp/crs_build_test_setup_20260224.log`
   - `/tmp/crs_check_ascran_test_setup_20260224.log` (`Status: 4 WARNINGs, 2 NOTEs`)

### 2026-02-24 - testthat deprecation/noise cleanup (`crsivderiv` tests)

Scope completed:

1. Removed deprecated `testthat::context()` usage from:
   - `/Users/jracine/Development/crs/tests/testthat/test-crsivderiv.R`
2. Added explicit warning suppression setting for deterministic test execution:
   - `display.warnings = FALSE` in `crsivderiv(...)` test calls.
3. Goal achieved:
   - direct `testthat::test_dir(...)` runs no longer emit the context deprecation warning and are less noisy.

Validation artifacts:

1. Direct tests:
   - `/tmp/crs_test_crsivderiv_modern_20260224.out` (`CRSIVDERIV_TESTS_OK`)
   - `/tmp/crs_test_fullsuite_post_crsivderiv_20260224.out` (`FULL_TESTS_OK`)
2. Deterministic install:
   - `/tmp/crs_install_test_crsivderiv_modern_20260224.log`
3. Tarball-first:
   - `/tmp/crs_build_test_crsivderiv_modern_20260224.log`
   - `/tmp/crs_check_ascran_test_crsivderiv_modern_20260224.log` (`Status: 4 WARNINGs, 3 NOTEs`)

### 2026-02-24 - B/R1.3 native matrix-kernel regression tests

Scope completed:

1. Added matrix-kernel regression tests in:
   - `/Users/jracine/Development/crs/tests/testthat/test-native-matrix-kernels.R`
2. Coverage added for:
   - `tensor.prod.model.matrix` (`mgcv_tmm`) parity against manual row-wise Kronecker construction.
   - `glp.model.matrix` (`glp_model_tmm`) structural invariants and column-membership within the tensor-product basis.
   - single-input edge case (`list(A)`) parity for both matrix constructors.
3. Scope is test-only; no runtime/package API changes.

Validation artifacts:

1. Direct tests:
   - `/tmp/crs_test_native_matrix_kernels_20260224.out` (`NATIVE_KERNEL_TESTS_OK`)
   - `/tmp/crs_test_fullsuite_post_matrixkernels_20260224.out` (`FULL_TESTS_OK`)
2. Deterministic install + smoke:
   - `/tmp/crs_install_matrix_kernels_20260224.log`
   - `/tmp/crs_loop_matrix_kernels_smoke_20260224.out` (`LOOP_SMOKE_OK`)
3. Tarball-first:
   - `/tmp/crs_build_matrix_kernels_20260224.log`
   - `/tmp/crs_check_ascran_matrix_kernels_20260224.log` (`Status: 4 WARNINGs, 3 NOTEs`)

### 2026-02-24 - C/R3.1 partial `.C` -> `.Call` migration (`uniquecombs`)

Scope completed:

1. Added new `.Call` wrapper in:
   - `/Users/jracine/Development/crs/src/RuniqueCombs.c`
   - `SEXP crs_uniquecombs_call(SEXP x)` allocates outputs in native code and preserves `index` mapping semantics.
2. Registered new native symbol in:
   - `/Users/jracine/Development/crs/src/crs_init.c`
3. Migrated R callsite in:
   - `/Users/jracine/Development/crs/R/mgcv.R`
   - `uniquecombs()` now uses `.Call("crs_uniquecombs_call", ..., PACKAGE="crs")`.
4. Result:
   - `.C(` callsites in `R/` reduced from `3` to `2`, with no user-facing API change.

Validation artifacts:

1. Deterministic install:
   - `/tmp/crs_install_uniquecombs_call_20260224.log`
2. Direct tests + suite:
   - `/tmp/crs_test_uniquecombs_call_native_20260224.out` (`NATIVE_TESTS_OK`)
   - `/tmp/crs_test_fullsuite_uniquecombs_call_20260224.out` (`FULL_TESTS_OK`)
3. Focused smoke:
   - `/tmp/crs_loop_uniquecombs_call_smoke_20260224.out` (`LOOP_SMOKE_OK`)
4. Tarball-first:
   - `/tmp/crs_build_uniquecombs_call_20260224.log`
   - `/tmp/crs_check_ascran_uniquecombs_call_20260224.log` (`Status: 4 WARNINGs, 3 NOTEs`)

### 2026-02-24 - C/R3.2 partial `.C` -> `.Call` migration (`gsl.bs`/`bs.des`)

Scope completed:

1. Added `.Call` wrappers in:
   - `/Users/jracine/Development/crs/src/gsl_bspline.c`
   - `crs_gsl_bspline_call(...)`
   - `crs_gsl_bspline_deriv_call(...)`
2. Registered wrappers in:
   - `/Users/jracine/Development/crs/src/crs_init.c`
3. Migrated R callsites in:
   - `/Users/jracine/Development/crs/R/gsl_bspline.R` (`bs.des()`)
   from `.C("gsl_bspline*")` to `.Call("crs_gsl_bspline*_call", ...)` while preserving matrix layout semantics.
4. Result:
   - `.C(` callsites in `R/` reduced from `2` to `0`; `.Call(` callsites in `R/` increased accordingly.

Validation artifacts:

1. Deterministic install:
   - `/tmp/crs_install_gsl_call_20260224.log`
2. Direct tests + suite:
   - `/tmp/crs_test_gsl_call_native_20260224.out` (`TARGET_TESTS_OK`)
   - `/tmp/crs_test_fullsuite_gsl_call_20260224.out` (`FULL_TESTS_OK`)
3. Focused smoke:
   - `/tmp/crs_spline_gsl_call_smoke_20260224.out` (`SPLINE_SMOKE_OK`)
   - `/tmp/crs_npglp_gsl_call_smoke_20260224.out` (`NPGLP_SMOKE_OK`)
4. Tarball-first:
   - `/tmp/crs_build_gsl_call_20260224.log`
   - `/tmp/crs_check_ascran_gsl_call_20260224.log` (`Status: 4 WARNINGs, 3 NOTEs`)

### 2026-02-24 - A/R1.21 remove remaining vectorized `ifelse` usage

Scope completed:

1. Removed the final two `ifelse(...)` callsites in:
   - `/Users/jracine/Development/crs/R/util.R`
   - `/Users/jracine/Development/crs/R/spline.R`
2. Replaced with explicit vectorized logical indexing to preserve behavior:
   - sign-aware epsilon replacement in `NZD(...)`
   - SVD reciprocal-threshold mask in `svd_lm_fit(...)`.
3. Result:
   - `ifelse(` count in `R/` is now `0`.

Validation artifacts:

1. Parse + inventory:
   - `/tmp/crs_parse_no_ifelse_20260224.log` (`NO_IFELSE_PARSE_OK`)
   - inline inventory result: `ifelse 0`
2. Deterministic install and tests:
   - `/tmp/crs_install_no_ifelse_20260224.log`
   - `/tmp/crs_test_fullsuite_no_ifelse_20260224.out` (`FULL_TESTS_OK`)
3. Focused smoke:
   - `/tmp/crs_spline_no_ifelse_smoke_20260224.out` (`SPLINE_SMOKE_OK`)
   - `/tmp/crs_npglp_no_ifelse_smoke_20260224.out` (`NPGLP_SMOKE_OK`)
4. Tarball-first:
   - `/tmp/crs_build_no_ifelse_20260224.log`
   - `/tmp/crs_check_ascran_no_ifelse_20260224.log` (`Status: 4 WARNINGs, 3 NOTEs`)

### 2026-02-24 - C/R3.3 `.Call` parity checks against legacy `.C` (`bs.des`)

Scope completed:

1. Extended native-interface tests in:
   - `/Users/jracine/Development/crs/tests/testthat/test-native-interfaces.R`
2. Added a local legacy comparator helper (`legacy_bs_des_c`) that calls registered `.C` entry points:
   - `"gsl_bspline"`
   - `"gsl_bspline_deriv"`
3. Added parity assertions that new `.Call` path (`crs:::bs.des`) matches legacy `.C` outputs for:
   - uniform-knot basis (`deriv = 0`)
   - uniform-knot mixed derivatives (`deriv` varying by row)
   - quantile-knot basis (`deriv = 0`)
4. Result:
   - new `.Call` migration has direct regression coverage against legacy interface behavior while `.C` symbols remain available internally.

Validation artifacts:

1. Native parity tests:
   - `/tmp/crs_test_native_interfaces_parity_20260224.out` (`PASS 18, FAIL 0`)
2. Full test suite:
   - `/tmp/crs_test_fullsuite_parity_20260224.out` (`PASS 75, WARN 1, FAIL 0`)
3. Deterministic install:
   - `/tmp/crs_install_native_parity_20260224.log`
4. Tarball-first:
   - `/tmp/crs_build_native_parity_20260224.log`
   - `/tmp/crs_check_ascran_native_parity_20260224.log` (`Status: 4 WARNINGs, 3 NOTEs`)

### 2026-02-24 - A/R1.22 `do.call` and sequence cleanup in matrix-construction paths

Scope completed:

1. Updated explicit function dispatch and sequence construction in:
   - `/Users/jracine/Development/crs/R/glp.model.matrix.R`
   - replaced `do.call('order', lapply(ncol(z):1, ...))` with `do.call(order, lapply(rev(seq_len(ncol(z))), ...))`, plus a defensive `ncol(z) > 0L` guard.
2. Updated `expand.grid` dispatch in:
   - `/Users/jracine/Development/crs/R/matrix.combns.R`
   - now uses `do.call(base::expand.grid, ...)` explicitly.
3. Removed ambiguous third-argument `do.call` usage in:
   - `/Users/jracine/Development/crs/R/np.regression.glp.R`
   - now uses `do.call(base::expand.grid, c(degree.list, list(KEEP.OUT.ATTRS = FALSE)))`.

Validation artifacts:

1. Targeted tests:
   - `/tmp/crs_test_docall_targeted_20260224.out` (`PASS 23, WARN 1, FAIL 0`)
2. Full test suite:
   - `/tmp/crs_test_docall_full_20260224.out` (`PASS 75, WARN 1, FAIL 0`)
3. Deterministic install:
   - `/tmp/crs_install_docall_20260224.log`
4. Tarball-first:
   - `/tmp/crs_build_docall_20260224.log`
   - `/tmp/crs_check_docall_20260224.log` (`Status: 4 WARNINGs, 2 NOTEs`)

### 2026-02-24 - Forensic closure sweep (R + non-NOMAD C audit)

Scope completed:

1. Ran consolidated pattern sweep over `R/` and non-NOMAD `src/` to verify remaining modernization debt.
2. Confirmed zero live matches in `R/` for:
   - string `do.call("...")` or `do.call('...')`
   - `eval(parse(...))`
   - `.C(` callsites
   - `ifelse(...)`
   - `<<-`
3. Confirmed only comment-only `1:length(...)` matches remain (`R/clsd.R` comments), with no active code usage.
4. At closure-sweep time, remaining `eval(...)` callsites (`R/crs.R`, `R/util.R`, `R/stepCV.R`) were flagged for dedicated parity-harness replacement work.
5. Non-NOMAD `src/` risky-call scan found only:
   - commented debug `Rprintf` in `src/mgcv.c`
   - NOMAD-wrapper messaging in `src/snomadr.h` / `src/snomadr.cpp` (explicitly out-of-scope for this modernization exercise).

Validation artifacts:

1. Consolidated forensic audit log:
   - `/tmp/crs_forensic_sweep_20260224.log`

### 2026-02-24 - A/R2.5 remove remaining direct `eval(...)` in R layer

Scope completed:

1. Added parity-safe internal evaluator helpers in:
   - `/Users/jracine/Development/crs/R/util.R`
   - `.crs_resolve_call_head(...)`
   - `.crs_as_eval_env(...)`
   - `.crs_eval_call(...)`
2. Replaced remaining direct `eval(...)` callsites in:
   - `/Users/jracine/Development/crs/R/stepCV.R`
     - `addterm.default(...)` and `dropterm.default(...)` now evaluate `update(..., evaluate = FALSE)` call objects via `.crs_eval_call(...)`.
   - `/Users/jracine/Development/crs/R/crs.R`
     - `crs.sigtest` call-data resolution now uses `.crs_eval_call(...)`.
   - `/Users/jracine/Development/crs/R/util.R`
     - `succeedWithResponse(...)` now uses `.crs_eval_call(...)`.
3. Added dedicated regression coverage in:
   - `/Users/jracine/Development/crs/tests/testthat/test-eval-helper.R`
   - parity against `eval(...)` for `update()` call-object evaluation,
   - namespace-qualified callable resolution,
   - response-presence detection behavior in `succeedWithResponse(...)`.
4. Result:
   - `eval(` count in `/Users/jracine/Development/crs/R` is now `0`.

Validation artifacts:

1. Forced clean install to avoid stale namespace bytecode during direct-test workflow:
   - `/tmp/crs_install_preclean_eval_refactor_20260224.log`
2. Targeted tests:
   - `/tmp/crs_test_eval_targeted_20260224.out` (`PASS 28, WARN 1, FAIL 0`)
3. Full test suite:
   - `/tmp/crs_test_eval_full_20260224.out` (`PASS 80, WARN 1, FAIL 0`)
4. Deterministic install:
   - `/tmp/crs_install_eval_refactor_20260224.log`
5. Tarball-first:
   - `/tmp/crs_build_eval_refactor_20260224.log`
   - `/tmp/crs_check_eval_refactor_20260224.log` (`Status: 4 WARNINGs, 4 NOTEs`)
   - `/tmp/crs_build_eval_refactor_clean_20260224.log`
   - `/tmp/crs_check_eval_refactor_clean_20260224.log` (`Status: 4 WARNINGs, 4 NOTEs`)

### 2026-02-24 - A/R1.23 `vapply` numeric-column detection sweep

Scope completed:

1. Replaced repeated `sapply(seq_len(ncol(...)), is.numeric)` patterns with `vapply(..., is.numeric, logical(1L))` in:
   - `/Users/jracine/Development/crs/R/util.R`
   - `/Users/jracine/Development/crs/R/np.regression.glp.R`
2. Simplified derived counters/indexes to scalar-safe direct forms:
   - `which(xdat.numeric == TRUE)` -> `which(xdat.numeric)`
   - `sum(... == TRUE)` -> `sum(xdat.numeric)`
   - `which(xdat.numeric == FALSE)` -> `which(!xdat.numeric)`
3. Removed final active `1:NROW(...)` indexing form in:
   - `/Users/jracine/Development/crs/R/glp.model.matrix.R`
   - `rownames(z) <- seq_len(NROW(z))`.
4. Result:
   - no remaining active `1:length(...)` / `1:ncol(...)` / `1:NCOL(...)` / `1:NROW(...)` loop-header/index forms in `R/`.

Validation artifacts:

1. Deterministic preclean install:
   - `/tmp/crs_install_vapply_sweep_20260224.log`
2. Targeted tests:
   - `/tmp/crs_test_vapply_targeted_20260224.out` (`PASS 33, WARN 1, FAIL 0`)
3. Full test suite:
   - `/tmp/crs_test_vapply_full_20260224.out` (`PASS 80, WARN 1, FAIL 0`)
4. Tarball-first:
   - `/tmp/crs_build_vapply_sweep_20260224.log`
   - `/tmp/crs_check_vapply_sweep_20260224.log` (`Status: 4 WARNINGs, 3 NOTEs`)

### 2026-02-24 - A/R1.24 boolean control simplification (`npglpreg` paths)

Scope completed:

1. Simplified boolean control branches in:
   - `/Users/jracine/Development/crs/R/np.regression.glp.R`
2. Replaced scalar/vector boolean comparisons with direct logical forms:
   - `xdat.numeric[i] == TRUE` -> `xdat.numeric[i]`
   - `xdat.numeric[i] != TRUE` -> `!xdat.numeric[i]`
   - `xdat.unordered[i] == TRUE` -> `xdat.unordered[i]`
   - `which(xdat.numeric == TRUE)` -> `which(xdat.numeric)`
3. Hardened scalar flag checks:
   - `leave.one.out == TRUE` -> `leave.one.out`
   - `mean == TRUE` -> `isTRUE(mean)`
4. Extended `vapply` use for predictor-type detection in additional local paths:
   - `xdat.numeric <- vapply(xdat, is.numeric, logical(1L))`
   - `xdat.unordered <- vapply(xdat, function(col) is.factor(col) && !is.ordered(col), logical(1L))`

Validation artifacts:

1. Deterministic preclean install:
   - `/tmp/crs_install_bool_sweep_20260224.log`
2. Targeted tests:
   - `/tmp/crs_test_bool_targeted_20260224.out` (`PASS 28, WARN 1, FAIL 0`)
3. Full test suite:
   - `/tmp/crs_test_bool_full_20260224.out` (`PASS 80, WARN 1, FAIL 0`)
4. Tarball-first:
   - `/tmp/crs_build_bool_sweep_20260224.log`
   - `/tmp/crs_check_bool_sweep_20260224.log` (`Status: 4 WARNINGs, 4 NOTEs`)
   - `/tmp/crs_build_bool_sweep_clean_20260224.log`
   - `/tmp/crs_check_bool_sweep_clean_20260224.log` (`Status: 4 WARNINGs, 4 NOTEs`)

### 2026-02-24 - A/R1.25 scalar logical gate cleanup (`crs` front-end)

Scope completed:

1. Simplified remaining scalar `== TRUE` checks in:
   - `/Users/jracine/Development/crs/R/crs.R`
2. Changes:
   - `kernel == TRUE` -> `isTRUE(kernel)`
   - `prune == TRUE` -> `isTRUE(prune)`
3. Updated branches:
   - categorical-kernel fallback when no factors are present (`num.z` check).
   - prune incompatibility checks (`kernel && prune`, `tau && prune`).
4. Result:
   - no active `== TRUE` / `!= TRUE` / `== FALSE` / `!= FALSE` in executable `R/` code (remaining matches are comment text only).

Validation artifacts:

1. Deterministic preclean install:
   - `/tmp/crs_install_bool_kernel_prune_20260224.log`
2. Targeted tests:
   - `/tmp/crs_test_bool_kernel_prune_targeted_20260224.out` (`PASS 27, WARN 1, FAIL 0`)
3. Full test suite:
   - `/tmp/crs_test_bool_kernel_prune_full_20260224.out` (`PASS 80, WARN 1, FAIL 0`)
4. Tarball-first:
   - `/tmp/crs_build_bool_kernel_prune_20260224.log`
   - `/tmp/crs_check_bool_kernel_prune_20260224.log` (`Status: 4 WARNINGs, 4 NOTEs`)

### 2026-02-24 - A/R1.26 centralized RNG seed-state helpers

Scope completed:

1. Added shared internal RNG seed helpers in:
   - `/Users/jracine/Development/crs/R/util.R`
   - `.crs_capture_seed(...)`
   - `.crs_restore_seed(...)`
2. Replaced duplicated seed save/restore blocks with helper calls in:
   - `/Users/jracine/Development/crs/R/clsd.R`
   - `/Users/jracine/Development/crs/R/snomadr.R`
   - `/Users/jracine/Development/crs/R/np.regression.glp.R`
   - `/Users/jracine/Development/crs/R/frscvNOMAD.R`
   - `/Users/jracine/Development/crs/R/krscvNOMAD.R`
   - `/Users/jracine/Development/crs/R/crssigtest.R`
3. Added helper regression tests in:
   - `/Users/jracine/Development/crs/tests/testthat/test-seed-helper.R`
4. Result:
   - no remaining `save.seed` / `exists.seed` duplication in `R/`.

Validation artifacts:

1. Deterministic install:
   - `/tmp/crs_install_seed_helper_20260224.log`
2. Targeted tests:
   - `/tmp/crs_test_seed_helper_targeted_20260224.out` (`PASS 13, WARN 1, FAIL 0`)
3. Full test suite:
   - `/tmp/crs_test_seed_helper_full_20260224.out` (`PASS 82, WARN 1, FAIL 0`)
4. Tarball-first:
   - `/tmp/crs_build_seed_helper_20260224.log`
   - `/tmp/crs_check_seed_helper_20260224.log` (`Status: 4 WARNINGs, 4 NOTEs`)

### 2026-02-24 - A/R1.27 scalar-loop index correctness in `npglpreg` paths

Scope completed:

1. Replaced scalar-loop headers that iterated only the last index in:
   - `/Users/jracine/Development/crs/R/np.regression.glp.R`
2. Key fixes:
   - `for(i in num.numeric)` -> `for(i in seq_len(num.numeric))`
   - `for(i in num.bw)` -> `for(i in seq_len(num.bw))`
3. Normalized adjacent degree-index slices to explicit index vectors:
   - `(num.bw+1):(num.bw+num.numeric)` -> `num.bw + seq_len(num.numeric)` in touched paths.
4. Behavioral intent:
   - preserve optimization behavior,
   - ensure warning/index checks are applied consistently across all predictors (not only the final one).

Validation artifacts:

1. Deterministic install:
   - `/tmp/crs_install_npglp_index_loops_20260224.log`
2. Targeted tests:
   - `/tmp/crs_test_npglp_index_loops_targeted_20260224.out` (`PASS 18, WARN 1, FAIL 0`)
3. Full test suite:
   - `/tmp/crs_test_npglp_index_loops_full_20260224.out` (`PASS 82, WARN 1, FAIL 0`)
4. Tarball-first:
   - `/tmp/crs_build_npglp_index_loops_20260224.log`
   - `/tmp/crs_check_npglp_index_loops_20260224.log` (`Status: 4 WARNINGs, 4 NOTEs`)

### 2026-02-24 - A/R1.28 registered-symbol `.Call` hygiene sweep

Scope completed:

1. Replaced remaining string-based `.Call("...")` invocations with registered-symbol calls (`.Call(symbol, ...)`) in:
   - `/Users/jracine/Development/crs/R/gsl_bspline.R`
   - `/Users/jracine/Development/crs/R/mgcv.R`
   - `/Users/jracine/Development/crs/R/spline.R`
   - `/Users/jracine/Development/crs/R/stepCV.R`
2. Removed now-redundant `PACKAGE = "crs"` callsite arguments on those paths.
3. Result:
   - `rg -n "\\.Call\\(\"" /Users/jracine/Development/crs/R` returns no matches.

Validation artifacts:

1. Deterministic install:
   - `/tmp/crs_install_callsymbols_20260224.log`
2. Targeted tests:
   - `/tmp/crs_test_callsymbols_targeted_20260224.out` (`PASS 36, WARN 0, FAIL 0`)
3. Full test suite:
   - `/tmp/crs_test_callsymbols_full_20260224.out` (`PASS 82, WARN 1, FAIL 0`)
4. Tarball-first:
   - `/tmp/crs_build_callsymbols_20260224.log`
   - `/tmp/crs_check_callsymbols_20260224.log` (`Status: 4 WARNINGs, 4 NOTEs`)

### 2026-02-24 - A/R1.29 robust `crs.messages` restoration on IV error paths

Scope completed:

1. Added function-scope option restoration guards in:
   - `/Users/jracine/Development/crs/R/crsiv.R`
   - `/Users/jracine/Development/crs/R/crsivderiv.R`
2. Change:
   - capture initial `crs.messages` at function entry,
   - restore with `on.exit(options(crs.messages = crs.messages), add = TRUE)` to prevent option leakage when errors occur while message-suppression toggles are active.
3. Added regression tests that force internal errors and assert option restoration:
   - `/Users/jracine/Development/crs/tests/testthat/test-crsiv.R`
   - `/Users/jracine/Development/crs/tests/testthat/test-crsivderiv.R`

Validation artifacts:

1. Deterministic install:
   - `/tmp/crs_install_option_restore_20260224.log`
2. Targeted tests:
   - `/tmp/crs_test_option_restore_targeted_20260224.out` (`PASS 20, WARN 0, FAIL 0`)
3. Full test suite:
   - `/tmp/crs_test_option_restore_full_20260224.out` (`PASS 86, WARN 1, FAIL 0`)
4. Tarball-first:
   - `/tmp/crs_build_option_restore_20260224.log`
   - `/tmp/crs_check_option_restore_20260224.log` (`Status: 4 WARNINGs, 4 NOTEs`)

### 2026-02-24 - A/R1.30 `crs.messages` unset-option guard in IV paths

Scope completed:

1. Hardened option capture/restore and scalar messaging gates in:
   - `/Users/jracine/Development/crs/R/crsiv.R`
   - `/Users/jracine/Development/crs/R/crsivderiv.R`
2. Changes:
   - store entry option value via `old.crs.messages <- getOption("crs.messages")`,
   - restore exact prior state on exit with `on.exit(options(crs.messages = old.crs.messages), add = TRUE)`,
   - use scalar-safe `crs.messages <- isTRUE(old.crs.messages)` for branch controls.
3. Added regression tests for unset option (`NULL`) behavior:
   - `/Users/jracine/Development/crs/tests/testthat/test-crsiv.R`
   - `/Users/jracine/Development/crs/tests/testthat/test-crsivderiv.R`
4. Result:
   - IV routines no longer risk `if (crs.messages)` errors when the option is unset,
   - option state is restored to its original value after both success and forced-error paths.

Validation artifacts:

1. Deterministic install:
   - `/tmp/crs_install_option_null_guard_20260224.log`
2. Targeted tests:
   - `/tmp/crs_test_option_null_guard_targeted_20260224.out` (`PASS 25, WARN 0, FAIL 0`)
3. Full test suite:
   - `/tmp/crs_test_option_null_guard_full_20260224.out` (`PASS 90, WARN 1, FAIL 0`)
4. Tarball-first:
   - `/tmp/crs_build_option_null_guard_20260224.log`
   - `/tmp/crs_check_option_null_guard_20260224.log` (`Status: 4 WARNINGs, 4 NOTEs`)

### 2026-02-24 - A/R1.31 eval-helper namespace-resolution hardening

Scope completed:

1. Hardened hidden namespace callable resolution in:
   - `/Users/jracine/Development/crs/R/util.R`
2. Change:
   - `.crs_resolve_call_head()` now uses `utils::getFromNamespace(...)` for `:::` heads, eliminating undefined global-function check noise while preserving behavior.
3. Expanded regression coverage in:
   - `/Users/jracine/Development/crs/tests/testthat/test-eval-helper.R`
4. Added explicit parity check for hidden namespace function calls:
   - `stats:::.MFclass(x)` route now exercised through `.crs_eval_call(...)`.
5. Result:
   - `R CMD check --as-cran` `R code for possible problems` gate is now `OK`,
   - overall check profile improved from `4 WARNINGs, 4 NOTEs` to `4 WARNINGs, 3 NOTEs`.

Validation artifacts:

1. Deterministic install:
   - `/tmp/crs_install_eval_nsqual_20260224.log`
2. Targeted tests:
   - `/tmp/crs_test_eval_nsqual_targeted_20260224.out` (`PASS 6, WARN 0, FAIL 0`)
3. Full test suite:
   - `/tmp/crs_test_eval_nsqual_full_20260224.out` (`PASS 91, WARN 1, FAIL 0`)
4. Tarball-first:
   - `/tmp/crs_build_eval_nsqual_20260224.log`
   - `/tmp/crs_check_eval_nsqual_20260224.log` (`Status: 4 WARNINGs, 3 NOTEs`)

### 2026-02-24 - A/R1.32 recursive source-artifact cleanup hardening

Scope completed:

1. Hardened package cleanup script for deep tree artifact removal:
   - `/Users/jracine/Development/crs/cleanup`
2. Replaced shallow glob deletion (`src/*.o`, `src/*/*.o`, etc.) with recursive file cleanup:
   - `find src -type f (...) -delete` for `.o`, `.so`, `.d`, `.dll`, `.a`, `.rc`,
   - `find src -type f -name 'Makedeps' -delete`.
3. Added explicit source-artifact ignore guards:
   - `/Users/jracine/Development/crs/.Rbuildignore`
   - `^src/.*\\.o$`
   - `^src/.*\\.so$`
4. Verified shell lint:
   - `shellcheck /Users/jracine/Development/crs/cleanup` (`clean`).
5. Result:
   - recursive cleanup removes deep NOMAD object trees (`225` files -> `0`),
   - source tarball no longer bundles `.o`/`.so` artifacts (`0` matches in `tar tzf`),
   - `R CMD check --as-cran` moved from object-file NOTE to clean source-package gate.

Validation artifacts:

1. Install/build/check path before recursive cleanup fix:
   - `/tmp/crs_install_rbuildignore_objs_20260224.log`
   - `/tmp/crs_build_rbuildignore_objs_20260224.log`
   - `/tmp/crs_check_rbuildignore_objs_20260224.log` (`Status: 4 WARNINGs, 3 NOTEs`)
2. Recursive cleanup verification:
   - local count check: `find src -type f \\( -name '*.o' -o -name '*.so' \\) | wc -l` (`225 -> 0` after `sh ./cleanup`)
3. Build/check after recursive cleanup:
   - `/tmp/crs_build_cleanup_recursive_20260224.log`
   - `/tmp/crs_check_cleanup_recursive_20260224.log` (`Status: 5 WARNINGs, 1 NOTE`)
4. Important note:
   - the additional warning is a pre-existing NOMAD compile warning now surfaced when objects are no longer bundled in source (`QPModelUtils.cpp:462`, bitwise/logical warning), not a modernization regression in `R/` or non-NOMAD C logic.

### 2026-02-24 - A/R1.33 `stepCV` loop-index hardening (empty-scope safety)

Scope completed:

1. Hardened loop headers in:
   - `/Users/jracine/Development/crs/R/stepCV.R`
2. Changes:
   - `for(i in seq(ns))` -> `for(i in seq_len(ns))` in `addterm.default` and `dropterm.default`,
   - `for(i in 1L:ns)` -> `for(i in seq_len(ns))` in `dropterm.glm`,
   - `setdiff(seq(ncol(x)), ii)` -> `setdiff(seq_len(ncol(x)), ii)` in `dropterm.glm`.
3. Added regression tests:
   - `/Users/jracine/Development/crs/tests/testthat/test-stepcv.R`
   - `dropterm.default(..., scope = character(0))` now executes without loop-index artifacts,
   - `dropterm.glm(..., scope = character(0))` now executes without loop-index artifacts.
4. Result:
   - fixed latent `seq(0)` / `1:0` traversal issues in empty-scope paths,
   - no behavioral regression on full test suite.

Validation artifacts:

1. Deterministic install:
   - `/tmp/crs_install_stepcv_seq_len_20260224.log`
2. Targeted tests:
   - `/tmp/crs_test_stepcv_dir_20260224.out` (`PASS 7, WARN 0, FAIL 0`)
3. Full test suite:
   - `/tmp/crs_test_stepcv_full_20260224.out` (`PASS 98, WARN 1, FAIL 0`)
4. Runtime smoke:
   - `/tmp/crs_stepcv_empty_scope_smoke_20260224.out` (`OK`)
5. Tarball-first:
   - `/tmp/crs_build_stepcv_seq_len_20260224.log`
   - `/tmp/crs_check_stepcv_seq_len_20260224.log` (`Status: 5 WARNINGs, 1 NOTE`)

### 2026-02-24 - A/R1.34 `snomadr` argument-validation hardening

Scope completed:

1. Hardened user-argument validation internals in:
   - `/Users/jracine/Development/crs/R/snomadr.R`
2. Changes:
   - replaced brittle `names(fargs)[2:length(fargs)]` with scalar-safe `names(fargs)[-1L]`,
   - simplified missing-argument failure path to report first missing required argument deterministically,
   - simplified extraneous-argument failure path to report first unexpected argument deterministically.
3. Added regression tests:
   - `/Users/jracine/Development/crs/tests/testthat/test-snomadr-args.R`
   - missing required `...` argument for `eval.f`,
   - extraneous `...` argument detection when `eval.f` formal requirements are supplied.
4. Result:
   - removed length/indexing footgun in argument-name slicing,
   - made user-facing argument-check errors deterministic and explicit.

Validation artifacts:

1. Deterministic install:
   - `/tmp/crs_install_snomadr_args_20260224.log`
2. Targeted tests:
   - `/tmp/crs_test_snomadr_args_targeted_20260224.out` (`PASS 2, WARN 0, FAIL 0`)
3. Full test suite:
   - `/tmp/crs_test_snomadr_args_full_20260224.out` (`PASS 100, WARN 1, FAIL 0`)
4. Tarball-first:
   - `/tmp/crs_build_snomadr_args_20260224.log`
   - `/tmp/crs_check_snomadr_args_20260224.log` (`Status: 5 WARNINGs, 1 NOTE`)

### 2026-02-24 - A/R1.35 stepCV call-evaluation helper alignment

Scope completed:

1. Replaced remaining `eval.parent(...)` path in:
   - `/Users/jracine/Development/crs/R/stepCV.R`
2. Changes:
   - `fit <- eval.parent(fit)` -> `fit <- .crs_eval_call(fit, parent.frame())` in the main `stepCV` update loop,
   - final model/keep list slicing now uses `seq_len(nm)` in place of `seq(nm)`.
3. Expanded regression coverage:
   - `/Users/jracine/Development/crs/tests/testthat/test-stepcv.R`
   - added `stepCV` update-path smoke test that exercises call-object evaluation and output structure.
4. Result:
   - removed the last active `eval.parent(...)` callsite from `R/`,
   - aligned step-model evaluation semantics with the centralized eval helper used elsewhere.

Validation artifacts:

1. Deterministic install:
   - `/tmp/crs_install_stepcv_evalcall_20260224.log`
2. Targeted tests:
   - `/tmp/crs_test_stepcv_evalcall_targeted_20260224.out` (`PASS 9, WARN 0, FAIL 0`)
3. Full test suite:
   - `/tmp/crs_test_stepcv_evalcall_full_20260224.out` (`PASS 102, WARN 1, FAIL 0`)
4. Tarball-first:
   - `/tmp/crs_build_stepcv_evalcall_20260224.log`
   - `/tmp/crs_check_stepcv_evalcall_20260224.log` (`Status: 5 WARNINGs, 0 NOTEs`)

### 2026-02-24 - A/R1.36 snomadr single-argument `...` guard hardening

Scope completed:

1. Hardened `snomadr` preflight argument checks for `eval.f(x)` signatures in:
   - `/Users/jracine/Development/crs/R/snomadr.R`
2. Changes:
   - normalized supplied-argument name handling (`NULL` names -> explicit placeholders),
   - added explicit guard: when `eval.f` has only `x`, any supplied `...` now fails immediately with deterministic argument-specific error.
3. Expanded regression coverage:
   - `/Users/jracine/Development/crs/tests/testthat/test-snomadr-args.R`
   - new test for extraneous argument rejection when `eval.f` has a single formal (`x`) only.
4. Result:
   - prevents repeated downstream optimizer-evaluation errors from undeclared extra arguments,
   - makes failure fast and user-facing at input-validation stage.

Validation artifacts:

1. Deterministic install:
   - `/tmp/crs_install_snomadr_singlearg_guard_20260224.log`
2. Targeted tests:
   - `/tmp/crs_test_snomadr_singlearg_guard_targeted_20260224.out` (`PASS 3, WARN 0, FAIL 0`)
3. Full test suite:
   - `/tmp/crs_test_snomadr_singlearg_guard_full_20260224.out` (`PASS 103, WARN 1, FAIL 0`)
4. Tarball-first:
   - `/tmp/crs_build_snomadr_singlearg_guard_20260224.log`
   - `/tmp/crs_check_snomadr_singlearg_guard_20260224.log` (`Status: 5 WARNINGs, 1 NOTE`)

### 2026-02-24 - A/R1.37 snomadr formal-matching parity for defaults and `...`

Scope completed:

1. Replaced legacy post-`x` argument matching in:
   - `/Users/jracine/Development/crs/R/snomadr.R`
2. Changes:
   - honors defaulted formals as optional (no longer treated as required),
   - supports functions declaring `...` by permitting additional named/unnamed arguments,
   - preserves strict rejection for extraneous args when `...` is absent,
   - preserves deterministic first-failure reporting for missing required formals and unexpected args.
3. Expanded regression coverage in:
   - `/Users/jracine/Development/crs/tests/testthat/test-snomadr-args.R`
   - new guard that defaulted arguments are not falsely flagged as required,
   - new guard that `...` routes accept extra args while still enforcing required non-default formals.
4. Result:
   - fixes false-positive "requires argument" failures for valid `eval.f` signatures with defaults,
   - aligns `snomadr` preflight checks with standard R formal semantics without touching NOMAD core paths.

Validation artifacts:

1. Test summary (targeted + full):
   - `/tmp/crs_test_snomadr_argmatch_summary_20260224.txt`
   - targeted: `PASS 5, WARN 0, FAIL 0`
   - full: `PASS 105, WARN 1, FAIL 0` (one pre-existing `npglpreg` warning)
2. Tarball-first:
   - `/tmp/crs_build_snomadr_argmatch_20260224.log`
   - `/tmp/crs_check_snomadr_argmatch_20260224.log` (`Status: 5 WARNINGs, 1 NOTE`)

### 2026-02-24 - A/R1.38 forensic legacy-marker sweep (R + non-NOMAD C)

Scope completed:

1. Ran a fresh static forensic sweep across `R/` and non-NOMAD `src/` files to verify no active legacy constructs remain.
2. Confirmed active R-layer modernization invariants remain intact:
   - `eval(parse(...))`, `eval(...)`, `eval.parent(...)`, `ifelse(...)`, string `do.call("...")`, `<<-`, and R-layer `.C(` callsites are all `0`.
3. Cleaned residual audit-marker comments in touched non-NOMAD surfaces:
   - `/Users/jracine/Development/crs/src/RuniqueCombs.c`
   - `/Users/jracine/Development/crs/R/np.regression.glp.R`
4. Result:
   - no active legacy markers remain in R or non-NOMAD C targets from this modernization scope,
   - remaining `1:length(` hits are comment-only historical notes (no executable callsites).

Validation artifacts:

1. Static sweep report:
   - `/tmp/crs_forensic_static_sweep_20260224.txt`
2. Syntax gate for touched R file:
   - `/tmp/crs_parse_npglpreg_comment_cleanup_20260224.out` (`PARSE_OK`)

### 2026-02-24 - A/R1.39 residual scalar bitwise-operator cleanup

Scope completed:

1. Replaced remaining scalar bitwise condition operators in control flow with scalar logical operators in:
   - `/Users/jracine/Development/crs/R/clsd.R`
   - `/Users/jracine/Development/crs/R/crs.R`
   - `/Users/jracine/Development/crs/R/stepCV.R`
   - `/Users/jracine/Development/crs/R/np.regression.glp.R`
2. Changes:
   - scalar `|` -> `||` in range-warning checks and direction flags,
   - scalar `&` -> `&&` in nullable-input guards and scalar degree checks.
3. Result:
   - removes remaining scalar bitwise footguns in actively executed control paths,
   - preserves vectorized `|`/`&` uses in matrix/index computations where element-wise logic is required.

Validation artifacts:

1. Syntax gate:
   - `/tmp/crs_parse_scalar_logic_cleanup_20260224.out` (`PARSE_OK`)
2. Targeted tests:
   - `/tmp/crs_test_scalar_logic_targeted_20260224.out` (`PASS 30, WARN 1, FAIL 0`; warning pre-existing in `test-npglpreg.R`)
3. Full test suite:
   - `/tmp/crs_test_scalar_logic_full_20260224.out` (`PASS 105, WARN 1, FAIL 0`)
4. Tarball-first:
   - `/tmp/crs_build_scalar_logic_cleanup_20260224.log`
   - `/tmp/crs_check_scalar_logic_cleanup_20260224.log` (`Status: 5 WARNINGs, 1 NOTE`)

### 2026-02-24 - A/R1.40 `plot.clsd` scalar logical follow-up

Scope completed:

1. Closed remaining scalar bitwise operators in `plot.clsd`:
   - `/Users/jracine/Development/crs/R/clsd.R`
2. Changes:
   - `!distribution & !derivative` -> `!distribution && !derivative` in both `er` branches.
3. Result:
   - removes final scalar bitwise control in the plotting path while preserving branch behavior.

Validation artifacts:

1. Syntax gate:
   - `/tmp/crs_parse_clsd_scalar_logic_followup_20260224.out` (`PARSE_OK`)
2. Full test suite:
   - `/tmp/crs_test_clsd_scalar_logic_followup_full_20260224.out` (`PASS 105, WARN 1, FAIL 0`)
3. Tarball-first:
   - `/tmp/crs_build_clsd_scalar_logic_followup_20260224.log`
   - `/tmp/crs_check_clsd_scalar_logic_followup_20260224.log` (`Status: 5 WARNINGs, 1 NOTE`)

### 2026-02-24 - A/R1.41 `W.glp` invariant-hoist micro-optimization

Scope completed:

1. Hoisted invariant degree maximum computation in `W.glp`:
   - `/Users/jracine/Development/crs/R/np.regression.glp.R`
2. Changes:
   - compute `degree.max <- max(degree)` once,
   - reuse `degree.max` for index filtering and loop guards.
3. Result:
   - removes repeated invariant recomputation inside polynomial-index setup,
   - no behavioral change (indexing logic and filtering remain identical).

Validation artifacts:

1. Syntax gate:
   - `/tmp/crs_parse_npglpreg_degree_max_hoist_20260224.out` (`PARSE_OK`)
2. Full test suite:
   - `/tmp/crs_test_npglpreg_degree_max_hoist_full_20260224.out` (`PASS 105, WARN 1, FAIL 0`)
3. Tarball-first:
   - `/tmp/crs_build_npglpreg_degree_max_hoist_20260224.log`
   - `/tmp/crs_check_npglpreg_degree_max_hoist_final_20260224.log` (`Status: 5 WARNINGs, 1 NOTE`)
4. Check rerun notes:
   - `/tmp/crs_check_npglpreg_degree_max_hoist_20260224.log` ended with OS-level kill during install (`Killed: 9`), matching the known sporadic pre-existing check instability.
   - `/tmp/crs_check_npglpreg_degree_max_hoist_rerun_20260224.log` failed examples due missing installed `np` in this session; resolved by local install from `/Users/jracine/Development/np-master` (`/tmp/np_install_for_crs_check_20260224.log`) before the final successful check.

### 2026-02-24 - A/R1.42 explicit boolean constant normalization

Scope completed:

1. Normalized legacy short boolean constants in active argument positions:
   - `/Users/jracine/Development/crs/R/util.R`
   - `/Users/jracine/Development/crs/R/np.regression.glp.R`
   - `/Users/jracine/Development/crs/R/crssigtest.R`
   - `/Users/jracine/Development/crs/R/frscvNOMAD.R`
   - `/Users/jracine/Development/crs/R/krscvNOMAD.R`
2. Changes:
   - `lower.tail=F` -> `lower.tail=FALSE`,
   - `replace=T` -> `replace=TRUE`,
   - aligned sampled-argument comments to `replace=TRUE`.
3. Result:
   - removes remaining ambiguous boolean shorthand in runtime callsites,
   - improves readability and static-audit clarity with no behavioral change.

Validation artifacts:

1. Syntax gate:
   - `/tmp/crs_parse_bool_constants_cleanup_20260224.out` (`PARSE_OK`)
2. Targeted tests:
   - `/tmp/crs_test_bool_constants_targeted_20260224.out` (`PASS 19, WARN 1, FAIL 0`)
3. Full test suite:
   - `/tmp/crs_test_bool_constants_full_20260224.out` (`PASS 105, WARN 1, FAIL 0`)
4. Tarball-first:
   - `/tmp/crs_build_bool_constants_cleanup_20260224.log`
   - `/tmp/crs_check_bool_constants_cleanup_20260224.log` (`Status: 5 WARNINGs`)
5. Check environment note:
   - first build attempt failed because `np` was unavailable in library (`ERROR: dependency 'np' is not available`);
   - resolved by local install from `/Users/jracine/Development/np-master` (`/tmp/np_install_for_crs_build_20260224.log`) before final successful build/check.

### 2026-02-24 - A/R1.43 duplicate helper-definition consolidation

Scope completed:

1. Removed redundant top-level helper definitions while preserving runtime behavior:
   - deleted duplicate `scale_robust` definition from:
     - `/Users/jracine/Development/crs/R/np.regression.glp.R`
   - deleted dead shadowed `is.fullrank` definition (early file copy) from:
     - `/Users/jracine/Development/crs/R/util.R`
2. Runtime behavior preserved by retaining the effective definitions:
   - `scale_robust` remains in `/Users/jracine/Development/crs/R/util.R`,
   - `is.fullrank` remains in `/Users/jracine/Development/crs/R/util.R` (the implementation that was already effective at runtime due later redefinition order).
3. Result:
   - removes duplicate-maintenance drift risk for shared numerical helpers,
   - clarifies a single source of truth for each helper without touching NOMAD core paths.

Validation artifacts:

1. Syntax gate:
   - `/tmp/crs_parse_helper_dedup_clean_20260224.out` (`PARSE_OK`)
2. Deterministic install before installed-namespace tests:
   - `/tmp/crs_install_helper_dedup_20260224.log`
3. Targeted tests:
   - `/tmp/crs_test_helper_dedup_targeted_20260224.out` (`PASS 31, WARN 1, FAIL 0`)
4. Full test suite:
   - `/tmp/crs_test_helper_dedup_full_postinstall_20260224.out` (`PASS 105, WARN 1, FAIL 0`)
5. Tarball-first:
   - `/tmp/crs_build_helper_dedup_20260224.log`
   - `/tmp/crs_check_helper_dedup_20260224.log` (`Status: 5 WARNINGs, 1 NOTE`)

### 2026-02-24 - A/R1.44 `npglpreg` option-restore stacking safety

Scope completed:

1. Hardened option restoration in:
   - `/Users/jracine/Development/crs/R/np.regression.glp.R`
2. Changes:
   - `on.exit(options(crs.messages = old.crs.messages))` ->
     `on.exit(options(crs.messages = old.crs.messages), add = TRUE)`.
3. Result:
   - preserves current restore behavior while ensuring composition with any additional `on.exit` handlers in the same execution path,
   - reduces accidental handler clobbering risk in future refactors/extensions.

Validation artifacts:

1. Syntax gate:
   - `/tmp/crs_parse_onexit_add_true_20260224.out` (`PARSE_OK`)
2. Deterministic install:
   - `/tmp/crs_install_onexit_add_true_20260224.log`
3. Targeted tests:
   - `/tmp/crs_test_onexit_add_true_targeted_20260224.out` (`PASS 13, WARN 1, FAIL 0`)
4. Full test suite:
   - `/tmp/crs_test_onexit_add_true_full_20260224.out` (`PASS 105, WARN 1, FAIL 0`)
5. Tarball-first:
   - `/tmp/crs_build_onexit_add_true_20260224.log`
   - `/tmp/crs_check_onexit_add_true_20260224.log` (`Status: 5 WARNINGs, 1 NOTE`)

### 2026-02-24 - A/R1.45 GLP sequence-construction edge hardening

Scope completed:

1. Hardened degenerate-dimension sequence creation in:
   - `/Users/jracine/Development/crs/R/glp.model.matrix.R`
2. Changes:
   - replaced `sets <- 1:dimen[1]` with `sets <- seq_len(dimen[1])` in `construct.tensor.prod`.
3. Result:
   - avoids malformed `1:0`-style range creation when the leading dimension is zero,
   - preserves behavior for all positive dimensions.

Validation artifacts:

1. Syntax gate:
   - `/tmp/crs_parse_glp_seq_len_20260224.out` (`PARSE_OK`)
2. Deterministic install:
   - `/tmp/crs_install_glp_seq_len_20260224.log`
3. Targeted tests:
   - `/tmp/crs_test_glp_seq_len_targeted_20260224.out` (`PASS 23, WARN 1, FAIL 0`)
4. Full test suite:
   - `/tmp/crs_test_glp_seq_len_full_20260224.out` (`PASS 105, WARN 1, FAIL 0`)
5. Tarball-first:
   - `/tmp/crs_build_glp_seq_len_20260224.log`
   - `/tmp/crs_check_glp_seq_len_20260224.log` (`Status: 5 WARNINGs, 1 NOTE`)

### 2026-02-24 - A/R1.46 GLP zero-dimension-safe index/rep hardening

Scope completed:

1. Hardened additional sequence/index construction in:
   - `/Users/jracine/Development/crs/R/glp.model.matrix.R`
2. Changes:
   - replaced residual `1:(...)` / `1:d2` patterns with `seq_len(...)`,
   - guarded first-block row initialization on `d2 > 0L`,
   - made replicated level construction explicit via `seq.int(0L, d2-1L)` and `seq.int(0L, i-1L)`,
   - normalized `for (i in 2:k)` loops to `seq.int(2L, k)` under existing `k > 1` guard.
3. Result:
   - removes remaining `1:0`-style edge-case hazards in GLP index generation,
   - keeps behavior identical for positive dimensions used in normal paths.

Validation artifacts:

1. Syntax gate:
   - `/tmp/crs_parse_glp_seq_hardening_20260224.out` (`PARSE_OK`)
2. Deterministic install:
   - `/tmp/crs_install_glp_seq_hardening_20260224.log`
3. Targeted tests:
   - `/tmp/crs_test_glp_seq_hardening_targeted_20260224.out` (`PASS 23, WARN 1, FAIL 0`)
4. Full test suite:
   - `/tmp/crs_test_glp_seq_hardening_full_20260224.out` (`PASS 105, WARN 1, FAIL 0`)
5. Tarball-first:
   - `/tmp/crs_build_glp_seq_hardening_20260224.log`
   - `/tmp/crs_check_glp_seq_hardening_20260224.log` (`Status: 5 WARNINGs, 1 NOTE`)

### 2026-02-24 - A/R1.47 residual GLP/NPGLP `1:n` slice hardening

Scope completed:

1. Replaced residual `1:n`-style index slices in active internals:
   - `/Users/jracine/Development/crs/R/glp.model.matrix.R`
   - `/Users/jracine/Development/crs/R/np.regression.glp.R`
2. Changes:
   - matrix column slices `1:nc` -> `seq_len(nc)` in GLP index assembly,
   - bandwidth/solution extraction `x[1:num.bw]` -> `x[seq_len(num.bw)]` in NPGLP CV/AIC objective and final bandwidth extraction paths.
3. Result:
   - removes remaining `1:0` edge-case exposure in these helper slices,
   - preserves behavior for standard positive-dimension cases.

Validation artifacts:

1. Syntax gate:
   - `/tmp/crs_parse_glp_bw_seq_len_followup_20260224.out` (`PARSE_OK`)
2. Deterministic install:
   - `/tmp/crs_install_glp_bw_seq_len_followup_20260224.log`
3. Targeted tests:
   - `/tmp/crs_test_glp_bw_seq_len_followup_targeted_20260224.out` (`PASS 23, WARN 1, FAIL 0`)
4. Full test suite:
   - `/tmp/crs_test_glp_bw_seq_len_followup_full_20260224.out` (`PASS 105, WARN 1, FAIL 0`)
5. Tarball-first:
   - `/tmp/crs_build_glp_bw_seq_len_followup_20260224.log`
   - `/tmp/crs_check_glp_bw_seq_len_followup_20260224.log` (`Status: 5 WARNINGs, 1 NOTE`)

### 2026-02-24 - A/R1.48 `frscv/krscv` wrapper index-slice normalization

Scope completed:

1. Normalized remaining direct `1:num.*` index slices in wrapper paths:
   - `/Users/jracine/Development/crs/R/frscv.R`
   - `/Users/jracine/Development/crs/R/frscvNOMAD.R`
   - `/Users/jracine/Development/crs/R/krscv.R`
   - `/Users/jracine/Development/crs/R/krscvNOMAD.R`
2. Changes:
   - replaced direct slices such as `x[1:num.x]` with `x[seq_len(num.x)]`,
   - applied same normalization to upper-bound/start-value index slices and extracted solution vectors in NOMAD wrappers.
3. Result:
   - removes residual `1:0`-style edge-case risk in these wrapper indexing paths,
   - retains existing behavior for normal positive-dimension cases.

Validation artifacts:

1. Syntax gate:
   - `/tmp/crs_parse_fr_kr_seq_len_sweep_20260224.out` (`PARSE_OK`)
2. Deterministic install:
   - `/tmp/crs_install_fr_kr_seq_len_sweep_20260224.log`
3. Targeted tests:
   - `/tmp/crs_test_fr_kr_seq_len_sweep_targeted_20260224.out` (`PASS 17, WARN 1, FAIL 0`)
4. Full test suite:
   - `/tmp/crs_test_fr_kr_seq_len_sweep_full_20260224.out` (`PASS 105, WARN 1, FAIL 0`)
5. Tarball-first:
   - `/tmp/crs_build_fr_kr_seq_len_sweep_20260224.log`
   - `/tmp/crs_check_fr_kr_seq_len_sweep_20260224.log` (`Status: 5 WARNINGs, 1 NOTE`)

### 2026-02-24 - A/R1.49 final executable `1:n` index-slice sweep

Scope completed:

1. Closed remaining executable `1:num.*` / `[1:n]` index slices in:
   - `/Users/jracine/Development/crs/R/frscv.R`
   - `/Users/jracine/Development/crs/R/krscv.R`
   - `/Users/jracine/Development/crs/R/spline.R`
2. Changes:
   - normalized remaining matrix/vector slices to `seq_len(...)`,
   - replaced `htt_aug[1:n]` and `htt_all[1:n][zz]` with `seq_len(n)` forms.
3. Result:
   - static sweep now shows no active executable `1:num.*`, `1:length`, or `1:ncol` patterns in `R/` (remaining hits are comments only),
   - behavior unchanged under existing tests/checks.

Validation artifacts:

1. Syntax gate:
   - `/tmp/crs_parse_final_index_slice_sweep_20260224.out` (`PARSE_OK`)
2. Deterministic install:
   - `/tmp/crs_install_final_index_slice_sweep_20260224.log`
3. Targeted tests:
   - `/tmp/crs_test_final_index_slice_sweep_targeted_20260224.out` (`PASS 27, WARN 1, FAIL 0`)
4. Full test suite:
   - `/tmp/crs_test_final_index_slice_sweep_full_20260224.out` (`PASS 105, WARN 1, FAIL 0`)
5. Tarball-first:
   - `/tmp/crs_build_final_index_slice_sweep_20260224.log`
   - `/tmp/crs_check_final_index_slice_sweep_20260224.log` (`Status: 5 WARNINGs, 1 NOTE`)

### 2026-02-24 - A/R1.50 final `1:num.*` closure pass (wrappers + static audit)

Scope completed:

1. Closed the last executable `1:num.*` slice sites in:
   - `/Users/jracine/Development/crs/R/frscv.R`
   - `/Users/jracine/Development/crs/R/krscv.R`
   - `/Users/jracine/Development/crs/R/spline.R`
2. Changes:
   - normalized remaining wrapper/matrix slices to `seq_len(num.x)`,
   - retained existing bounds/shape semantics and explicit `drop=FALSE` guards where already present.
3. Result:
   - static pattern scan now reports no executable `1:num.*` or related legacy slice patterns in `R/` (remaining hits are comment-only in `clsd.R`),
   - no regressions in tests or tarball checks.

Validation artifacts:

1. Syntax gate:
   - `/tmp/crs_parse_final_seq_sweep2_20260224.out` (`PARSE_OK`)
2. Deterministic install:
   - `/tmp/crs_install_final_seq_sweep2_20260224.log`
3. Targeted tests:
   - `/tmp/crs_test_final_seq_sweep2_targeted_20260224.out` (`PASS 27, WARN 1, FAIL 0`)
4. Full test suite:
   - `/tmp/crs_test_final_seq_sweep2_full_20260224.out` (`PASS 105, WARN 1, FAIL 0`)
5. Tarball-first:
   - `/tmp/crs_build_final_seq_sweep2_20260224.log`
   - `/tmp/crs_check_final_seq_sweep2_20260224.log` (`Status: 5 WARNINGs, 1 NOTE`)

### 2026-02-24 - A/R1.51 guarded loop-header `2:n` normalization sweep

Scope completed:

1. Normalized remaining parity-safe guarded loop headers in:
   - `/Users/jracine/Development/crs/R/frscv.R`
   - `/Users/jracine/Development/crs/R/krscv.R`
   - `/Users/jracine/Development/crs/R/np.regression.glp.R`
   - `/Users/jracine/Development/crs/R/spline.R`
   - `/Users/jracine/Development/crs/R/util.R`
2. Changes:
   - replaced guarded `for (... in 2:n)` patterns with `for (... in seq.int(2L, n))` under existing `n > 1` checks,
   - preserved existing loop-body logic and ordering.
3. Result:
   - removes residual `2:1` range hazards in the remaining guarded active loops,
   - keeps behavior unchanged for positive-dimension paths already covered by current tests.

Validation artifacts:

1. Syntax gate:
   - `/tmp/crs_parse_loop_seqint_guarded_20260224.out` (`PARSE_OK`)
2. Deterministic install:
   - `/tmp/crs_install_loop_seqint_guarded_20260224.log`
3. Targeted tests:
   - `/tmp/crs_test_loop_seqint_guarded_targeted_20260224.out` (`PASS 41, WARN 1, FAIL 0`)
4. Full test suite:
   - `/tmp/crs_test_loop_seqint_guarded_full_20260224.out` (`PASS 105, WARN 1, FAIL 0`)
5. Tarball-first:
   - `/tmp/crs_build_loop_seqint_guarded_20260224.log`
   - `/tmp/crs_check_loop_seqint_guarded_20260224.log` (`Status: 5 WARNINGs, 1 NOTE`)

### 2026-02-24 - A/R1.52 IV iteration-loop guard hardening (`iterate.max`)

Scope completed:

1. Hardened residual unguarded `2:iterate.max` loop headers in:
   - `/Users/jracine/Development/crs/R/crsiv.R`
   - `/Users/jracine/Development/crs/R/crsivderiv.R`
2. Changes:
   - normalized loop headers to `if (iterate.max > 1L) for (j in seq.int(2L, iterate.max))`,
   - initialized `convergence <- "ITERATE_MAX"` before loop entry for deterministic no-loop edge behavior.
3. Result:
   - removes remaining `2:1` range hazards in active IV iterative paths,
   - preserves existing behavior for standard `iterate.max >= 2` use while making edge semantics explicit.

Validation artifacts:

1. Syntax gate:
   - `/tmp/crs_parse_iterate_seq_guard_clean_20260224.out` (`PARSE_OK`)
2. Deterministic install:
   - `/tmp/crs_install_iterate_seq_guard_20260224.log`
3. Targeted tests:
   - `/tmp/crs_test_iterate_seq_guard_targeted_20260224.out` (`PASS 28, WARN 0, FAIL 0`)
4. Full test suite:
   - `/tmp/crs_test_iterate_seq_guard_full_20260224.out` (`PASS 105, WARN 1, FAIL 0`)
5. Tarball-first:
   - `/tmp/crs_build_iterate_seq_guard_20260224.log`
   - `/tmp/crs_check_iterate_seq_guard_20260224.log` (`Status: 5 WARNINGs`)
