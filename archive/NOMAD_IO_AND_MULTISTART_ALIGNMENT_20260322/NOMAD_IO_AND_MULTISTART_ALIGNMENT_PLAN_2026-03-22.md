# NOMAD I/O And Multistart Alignment Plan (`crs`)

Date: 2026-03-22
Status: Proposed plan only; no `crs` code edits made in this investigation.

## Goal

Bring `crs` NOMAD behavior into alignment with the restored `np` contract:

1. Non-NOMAD routes keep their current user-visible I/O.
2. NOMAD routes use one managed progress line that remains visible and clears on finish.
3. For NOMAD routes, that one line carries rich search state:
   - current multistart,
   - evaluation/iteration count,
   - current degree,
   - best degree so far,
   - current/best objective value.
4. `nmulti` behavior for polynomial-order search should follow the same explicit outer-restart pattern now used in `npglpreg` and in `np-master` LP/NOMAD search.

## Repo State At Planning Time

Current `crs` branch and recent history:

1. Branch: `master`
2. Recent NOMAD-related commits on this line:
   - `01defb1` `Honor nmulti in crs default selection path`
   - `d28aa48` `Fix NOMAD handoff cleanup references`
   - `967cb52` `Polish NOMAD progress handoff`
   - `0d3dff8` `Fix single-line progress handoff cleanup`
   - `b9f10a1` `Unify NOMAD progress line ownership in crs`
   - `6a6b907` `Make NOMAD nmulti use explicit restarts`

Uncommitted worktree state observed in the main checkout:

1. Modified:
   - `CRS_MODERNIZATION_PLAN.md`
   - `DESCRIPTION`
   - `NOMAD_STATUS_AND_GATES.md`
2. Deleted:
   - `NOMAD_FORWARD_GUIDE.md`
   - `NOMAD_INTEGRATION_MAP.md`
3. Untracked:
   - `archive/NOMAD_MD_CONSOLIDATION_20260322/`

Attached worktrees observed:

1. Main checkout:
   - `/Users/jracine/Development/crs`
2. Detached historical/check worktrees:
   - `/private/tmp/crs-ascran-20260322`
   - `/private/tmp/crs_prev_nomad_check.QHjG0Q`
   - several older detached `crs_*` worktrees marked prunable

None of the above non-plan changes were edited by this planning pass.

## Current Functional Status

### 1. `npglpreg`

`npglpreg` is partially modernized already.

What is already in place:

1. Explicit outer restart handling in `glpcvNOMAD()` rather than relying only on `snomadr(..., nmulti=...)`.
2. Degree-start helpers modeled on the newer `np` approach:
   - `.crs_glp_nomad_build_degree_starts()`
   - `.crs_glp_nomad_build_start_matrix()`
3. Stored restart artifacts:
   - `nomad.restart.results`
   - `nomad.restart.fval`
   - `nomad.best.restart`
   - `nomad.starts`

What is still missing:

1. `np`-style rich one-line NOMAD progress.
2. Best-so-far/current-degree/current-eval detail on that line.
3. A proper managed handoff for Powell-like refinement messaging if we add that path later.

Observed current live I/O for `npglpreg(..., cv="degree-bandwidth", nmulti=3, ...)`:

1. `Calling NOMAD ...`
2. `Calling NOMAD ... multistart 1/3`
3. repeated `fv = ...`

This is much coarser than the restored `np` NOMAD line.

### 2. `crs()`

`crs()` is not on the same multistart structure as `npglpreg`.

Current routing:

1. `crs()` still delegates NOMAD CV to the older helpers:
   - `frscvNOMAD()`
   - `krscvNOMAD()`
2. Those helpers still pass `nmulti` down into `snomadr(..., nmulti=...)`.
3. They do not yet use the newer explicit outer-restart matrix logic that `npglpreg` uses.

Observed current live I/O for `crs(..., cv="nomad", nmulti=3, ...)`:

1. old NOMAD display stream with duplicated lines under the capture harness,
2. degree/segment/CV detail from the older `DISPLAY_DEGREE` path,
3. not the newer single managed line used in `np`.

### 3. `crsiv` and `crsivderiv`

These currently expose only outer IV-stage progress, not nested NOMAD detail.

Why:

1. Both functions repeatedly call `.crs_set_messages(..., FALSE)` around nested `crs()` fits.
2. They replace nested search visibility with coarse stage labels such as:
   - `E[y|z]`
   - `E[y-phi(z)|w]`
   - `E[E[y-phi(z)|w]|z]`
3. So even if `crs()` NOMAD I/O improves, `crsiv` and `crsivderiv` will still suppress it unless the ownership/handoff model changes.

## Proposed Fix

Use the `np` restoration as the design pattern, but apply it in `crs` in two layers:

### Layer A: Shared NOMAD progress payload for `crs`

Introduce a shared NOMAD progress state/detail formatter in `crs` that can express:

1. `multistart i/n`
2. `eval k`
3. `deg (...)` when degree is searchable
4. `best (...)`
5. `fv=...`
6. optional `iteration ...` if the route has it

This should be managed by the existing `crs` single-line progress system in `R/progress.R`, not by standalone messages.

### Layer B: Shared explicit outer multistart for all degree-search NOMAD routes

Standardize degree-search multistart around explicit outer restarts rather than mixed behavior.

Target direction:

1. Keep `npglpreg` as the reference implementation inside `crs`.
2. Port the same explicit outer-restart structure to:
   - `frscvNOMAD()`
   - `krscvNOMAD()`
3. Reuse one centralized degree-start builder rather than maintaining separate ad hoc `x0` generation.

This would align `crs()` with the same conceptual model now used in `npglpreg` and restored in `np-master`.

### Layer C: IV/derivative handoff instead of suppression

For `crsiv` and `crsivderiv`:

1. keep the outer IV-stage progress design,
2. but when a nested `crs()` call enters a NOMAD route, let the nested NOMAD owner take over the visible line,
3. then resume the outer IV-stage line after the nested fit completes.

That is the same ownership problem solved in `np`: the outer wrapper should not mask the richer inner NOMAD stream.

## Recommended Tranche Plan

### Tranche 1: `npglpreg` I/O only

Scope:

1. No search-logic changes.
2. No objective changes.
3. No multistart-semantics changes.
4. Only replace current coarse `Calling NOMAD` / `fv = ...` progress with one managed rich line.

Files likely touched:

1. `R/np.regression.glp.R`
2. `R/progress.R`
3. targeted `testthat` progress contracts for `npglpreg`

Success criterion:

1. `npglpreg(... cv="degree-bandwidth", nmulti > 1 ...)` shows one rich NOMAD line with `multistart`, `eval`, `deg`, `best`, and `fv`.
2. Non-NOMAD `npglpreg` output remains unchanged.

### Tranche 2: `crs()` multistart alignment

Scope:

1. Migrate `frscvNOMAD()` and `krscvNOMAD()` from C-managed `nmulti` to explicit outer restarts at the R layer.
2. Centralize start generation so polynomial-order starts follow the same explicit pattern as `npglpreg`.

Files likely touched:

1. `R/frscvNOMAD.R`
2. `R/krscvNOMAD.R`
3. shared helper location, likely `R/util.R` or a new NOMAD helper section in an existing file

Success criterion:

1. `crs(..., cv="nomad", nmulti > 1)` uses explicit outer restarts.
2. chosen starts and best-restart bookkeeping are inspectable in returned objects or diagnostics.
3. numerical behavior remains stable apart from accepted restart-order differences.

### Tranche 3: `crs()` NOMAD one-line I/O

Scope:

1. Put `frscvNOMAD()` and `krscvNOMAD()` on the same managed-line NOMAD detail format as `npglpreg`.

Success criterion:

1. `crs(..., cv="nomad")` NOMAD routes show one managed line, not the old duplicated `DISPLAY_DEGREE` stream.
2. line includes degree/current/best/multistart/fv detail where applicable.

### Tranche 4: `crsiv` / `crsivderiv` handoff only

Scope:

1. Do not change IV math.
2. Do not change stopping rules.
3. Only change progress ownership so nested NOMAD fits remain visible.

Files likely touched:

1. `R/crsiv.R`
2. `R/crsivderiv.R`
3. possibly `R/progress.R`

Success criterion:

1. outer IV stage line remains readable,
2. nested NOMAD fits surface the same rich NOMAD line as `crs()` / `npglpreg`,
3. control returns cleanly to the IV-stage line after nested fit completion.

## Validation Plan

For every tranche:

1. keep changes I/O-only unless the tranche explicitly targets multistart semantics,
2. save pre/post artifacts under dated `/tmp` paths,
3. compare live output using the existing progress shadow harness,
4. run focused installed-package smokes for:
   - `npglpreg`
   - `crs`
   - `crsiv`
   - `crsivderiv`
5. confirm no regression in silent/noninteractive behavior when `display.nomad.progress = FALSE`

Critical route checks:

1. `npglpreg(y ~ x, cv="degree-bandwidth", nmulti=3, ...)`
2. `crs(y ~ x, basis="additive", cv="nomad", nmulti=3, ...)`
3. `crsiv(..., basis="additive", iterate.max=2, display.nomad.progress=TRUE, nmulti=...)`
4. `crsivderiv(..., basis="additive", iterate.max=2, display.nomad.progress=TRUE, nmulti=...)`

## Recommended Implementation Order

1. `npglpreg` I/O first.
2. `crs()` multistart alignment second.
3. `crs()` NOMAD I/O third.
4. `crsiv` / `crsivderiv` handoff last.

Reason:

1. `npglpreg` already has the closest structure to the target.
2. `crs()` still has an older multistart core, so I/O and semantics are entangled there.
3. `crsiv` and `crsivderiv` should not be touched until the nested `crs()` NOMAD owner is trustworthy.

## Decision Summary

The safest path is not one giant rewrite.

The right fix is:

1. treat `npglpreg` as the near-term reference path inside `crs`,
2. migrate `crs()` off the older `snomadr(..., nmulti=...)` multistart pattern,
3. then let `crsiv` / `crsivderiv` hand off to the nested NOMAD owner rather than suppressing it.

This keeps the work consistent with the `np` restoration while respecting the fact that `crs` currently has two different NOMAD architectures in the same repo.
