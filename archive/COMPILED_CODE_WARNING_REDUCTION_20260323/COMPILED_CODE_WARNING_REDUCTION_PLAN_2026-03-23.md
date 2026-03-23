# Compiled-Code Warning Reduction Plan

Date: 2026-03-23
Repo: `/Users/jracine/Development/crs`
Scope: Reduce the `R CMD check --as-cran` compiled-code warning from vendored NOMAD/SGTELIB with the highest possible safety and the lowest possible risk of estimator, optimizer, or API regression.

## Goal

Reduce or eliminate the `compiled code` warning emitted by `R CMD check --as-cran` for [`crs/libs/crs.so`](/Users/jracine/Development/crs/libs/crs.so) while preserving:

- all existing estimator behavior,
- all existing NOMAD optimization behavior,
- all existing user-visible R-layer I/O behavior restored in the recent I/O tranche,
- all package APIs and defaults,
- all current numerical results up to floating-point noise only where unavoidable.

This is a release-candidate safety tranche, not a modernization tranche.

## Non-Negotiable Constraints

1. No estimator drift.
   - No change in selected models, objective values, stopping behavior, or returned objects attributable to this tranche.

2. No optimizer-semantic drift.
   - No change to NOMAD search strategy, SGTELIB activation semantics, `nmulti`, random-start logic, budgets, stopping rules, or parameter parsing unless explicitly authorized after evidence review.

3. No user-facing R I/O regression.
   - The recently repaired single-line and IV-progress behavior for `crs`, `npglpreg`, `crsiv`, and `crsivderiv` must remain intact.

4. No hidden feature removal.
   - Do not silently disable SGTELIB-backed runtime capability merely to reduce check warnings.

5. One risk axis only.
   - This tranche is about compiled-code warning reduction only. No opportunistic refactors, cleanup, or interface work.

## Problem Statement

`R CMD check --as-cran` currently reports compiled-code warnings because the package builds vendored NOMAD/SGTELIB sources that contain:

- `std::cout`
- `std::cerr`
- `printf` / `puts` / `putchar`
- `exit`
- `rand`

The current build wiring in [`src/Makevars`](/Users/jracine/Development/crs/src/Makevars) compiles most of the vendored NOMAD and SGTELIB source trees using `find ... '*.cpp'`, with `-DUSE_SGTELIB` enabled.

Important forensic conclusion:

- this is not only one dead test file,
- some warning sources are clearly non-runtime/debug baggage,
- but several flagged objects are genuine runtime-compiled sources,
- therefore a safe solution must be staged and evidence-driven.

## Risk Assessment

### Low-Risk Candidates

These are candidate changes that are plausibly behavior-neutral and should be considered first:

1. Excluding obviously non-runtime vendor sources from the build if they are not referenced by runtime package code.
   - Example candidate already identified: `ext/sgtelib/src/Tests.cpp`

2. Patching vendor debug or warning output sites that do not affect control flow or returned values.
   - Example pattern: replacing unconditional `std::cout`/`std::cerr` diagnostics with no-op or R-safe logging shim only when those lines are not part of runtime semantics.

### Medium-Risk Candidates

These may still be acceptable, but only after route proof:

1. Patching vendor exception-reporting `printf` calls in the C interface.
2. Replacing `rand()` with a deterministic local generator inside SGTELIB utility helpers if and only if that code path is proven runtime-relevant and behavior-neutral in package usage.

### High-Risk Candidates

These are not acceptable for the first checkpoint:

1. Disabling `-DUSE_SGTELIB` globally.
2. Removing whole SGTELIB model or quad-model compilation paths without feature proof.
3. Rewriting vendor logic related to search, model construction, QP solver behavior, or stopping.
4. Changing NOMAD C/C++ parameter semantics or callback flow.

## Staged Execution Plan

### Stage 0: Freeze Baseline

Before any source edit, record the exact current behavior and warning surface.

Artifacts to save:

- current `R CMD check --as-cran` log,
- current unresolved-symbol scan of `crs.so`,
- current `git status`,
- focused runtime smoke logs for NOMAD routes.

Baseline behavior to freeze:

1. `crs(...)` NOMAD route
2. `npglpreg(...)` NOMAD route
3. `crsiv(...)` route
4. `crsivderiv(...)` route

Baseline outputs to freeze:

- selected degree/segment choices,
- objective values where practical,
- single-line progress behavior,
- check warning inventory.

Acceptance gate:

- no edits until baseline artifacts exist and are readable.

### Stage 1: Dead-Weight Build Exclusion Audit

Audit vendor source files that are strong candidates for exclusion from package build without affecting runtime.

Initial candidates:

- `src/nomad4_src/ext/sgtelib/src/Tests.cpp`

For each candidate, prove:

1. it is not required by any runtime-linked package path,
2. excluding it does not break build,
3. excluding it does not change focused runtime outputs,
4. excluding it reduces at least part of the compiled-code warning.

Implementation rule:

- change only [`src/Makevars`](/Users/jracine/Development/crs/src/Makevars) in this stage.

Checkpoint decision:

- keep the exclusion only if build and behavior remain unchanged and warning surface is measurably reduced,
- otherwise revert immediately.

### Stage 2: Runtime Callsite Classification

For remaining flagged symbols after Stage 1, classify each source occurrence into one of:

1. diagnostic-only output,
2. exception/error reporting,
3. control-flow-affecting termination,
4. true algorithmic randomness,
5. unknown/risky.

No edits in this stage beyond comments or notes if absolutely necessary.

Required output of this stage:

- a short table mapping each surviving warning symbol to:
  - source file,
  - runtime relevance,
  - proposed safe action,
  - risk level.

Checkpoint rule:

- no runtime vendor patching until each surviving symbol is classified.

### Stage 3: Lowest-Risk Runtime Silencing

Patch only callsites classified as both:

- runtime-compiled, and
- semantically diagnostic-only.

Allowed examples:

- replacing bare `std::cout`/`std::cerr` diagnostics with suppressed output when they are not part of returned-state logic,
- removing debug-only text prints from vendor code paths.

Not allowed in this stage:

- changing thrown exception types,
- changing return codes,
- changing stop criteria,
- changing search/model-selection logic.

Checkpoint decision:

- if warning count drops and all gates hold, keep the patch,
- if any numerical, behavioral, or I/O drift appears, revert immediately.

### Stage 4: Termination and RNG Review

Address `_exit` / `rand` only if still necessary after earlier stages.

This stage is explicitly higher risk and must be treated as optional unless the earlier stages are insufficient.

For `_exit`:

- patch only if the path is reachable from package runtime,
- prefer replacing with exception propagation or ordinary error return only if runtime-equivalent behavior can be proven.

For `_rand`:

- patch only if the path is reachable in actual package execution,
- if unreachable, prefer source exclusion over semantic replacement,
- if reachable, replacement must be justified with strong parity evidence.

Checkpoint decision:

- do not proceed to this stage without reviewing earlier evidence first.

## Validation Gates

These are hard keep/revert gates for every stage.

### Build Gates

1. Package must still build from tarball.
2. `R CMD check --as-cran` must not gain any new warning, note, or error.
3. The compiled-code warning inventory must be reduced or unchanged only for a reverted failed attempt.

### Numerical/Behavior Gates

For touched checkpoints, rerun at least:

1. `crs(...)` NOMAD selector smoke
2. `npglpreg(...)` NOMAD selector smoke
3. `crsiv(...)` smoke
4. `crsivderiv(...)` smoke, to the extent the standing unrelated evaluation-data issue allows

Required invariants:

- no unexpected drift in selected degree/segment,
- no unexpected drift in objective values beyond floating-point noise,
- no changed warnings/errors/messages except those intentionally tied to compiled-code silencing.

### I/O Gates

1. `crs` and `npglpreg` must retain the rich one-line NOMAD progress behavior.
2. `crsiv` and `crsivderiv` must retain outer-IV progress ownership.
3. No restored vendor stdout/stderr chatter may leak into R output.

## Commit Strategy

Use small checkpoint commits only.

Planned sequence:

1. Stage 1 exclusion-only commit, if it proves safe and useful.
2. Stage 3 low-risk vendor-silencing commit, only if Stage 1 is insufficient and Stage 2 classification supports it.
3. Additional commit(s) only if a later stage is explicitly justified by evidence.

Each commit must be independently revertible and independently validated.

## Explicit No-Go Decisions

Do not do any of the following in the first pass:

- no bulk vendor rewrite,
- no migration away from vendored NOMAD,
- no disabling SGTELIB globally,
- no performance tuning,
- no Makevars cleanup unrelated to this warning,
- no R-layer feature work,
- no unrelated documentation churn.

## Decision Framework After First Checkpoint

After Stage 1 we should explicitly decide between:

1. Stop if the reduction is meaningful and safe enough for RC.
2. Continue into low-risk runtime silencing if remaining warnings are still substantial and the callsites are clearly diagnostic-only.
3. Defer `_exit` / `rand` work if it requires semantics-risky vendor patching.

## Deliverables

By the end of the first checkpoint, we should have:

1. a precise baseline artifact set,
2. a proof about whether `Tests.cpp` exclusion is safe,
3. a before/after compiled-code warning diff,
4. a recommendation on whether deeper vendor patching is justified.

That first checkpoint should be enough to decide whether the warning can be reduced safely near release candidate, or whether further work would carry disproportionate regression risk.
