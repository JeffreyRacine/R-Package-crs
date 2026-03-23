# Revised NOMAD I/O And Multistart Alignment Plan (`crs`)

Date: 2026-03-22
Status: Revised after engineering-risk review

## Further Tightening For Release-Candidate Safety

This section narrows the plan further to maximize success probability and minimize collateral damage.

### Non-Negotiable Constraints

1. Do not touch `src/snomadr.cpp`, `src/snomadr.h`, `src/nomad4_src/**`, or `R/snomadr.R` in any I/O tranche.
2. Do not change `nmulti`, `random.seed`, `x0`, restart ordering, or evaluation-budget semantics in any I/O tranche.
3. Do not introduce a new shared abstraction until one route-specific implementation has already succeeded and been validated.
4. Do not change more than one visible NOMAD route family in a single commit.
5. Do not mix progress-format changes with numerical or performance changes.

### Stronger Design Rule

Prefer route-local adapters first, shared helpers second.

Why:

1. premature commonization increases collateral-regression risk,
2. `npglpreg`, `crs()`, and `crsiv` do not currently share the same NOMAD path,
3. proving the design on one route before extraction is safer than guessing.

## Why The First Plan Was Too Risky

The earlier plan had the right destination but too much structural ambition for a release-candidate phase.

Main problems:

1. It mixed two risk axes too early:
   - NOMAD I/O restoration
   - multistart semantics unification
2. It treated `npglpreg` as the natural reference implementation inside `crs`, but `crs()` and the IV routes still depend on older NOMAD plumbing. Making `npglpreg` the internal reference too early could force structural churn instead of a surgical fix.
3. It proposed migrating `frscvNOMAD()` / `krscvNOMAD()` off `snomadr(..., nmulti=...)` before proving that this is required to restore the desired user-visible behavior.
4. It did not make a hard enough distinction between:
   - restoring user-visible NOMAD line ownership,
   - changing how restarts are generated and executed.
5. It risked violating the workspace rule of one risk axis per tranche.
6. It assumed the best fix was architectural convergence. At this point, the safest fix is behavioral restoration with minimum internal disturbance.

Under a highest-standards / lowest-regression lens, the right question is:

What is the smallest set of changes that can restore the exact user-visible NOMAD contract without changing solver semantics unless we prove they are broken?

## Revised Engineering Principles

These principles should govern the work:

1. Treat I/O restoration and multistart semantics as separate projects unless a hard blocker proves otherwise.
2. Do not replace working restart engines during an I/O tranche.
3. Prefer adapters and hooks over rewrites.
4. Add contract tests before edits whenever possible.
5. Use live output captures plus returned-object checks to prove no collateral change.
6. Do not standardize internal architecture just because two paths are similar; standardize only when the current divergence itself causes the bug.

## Revised Diagnosis

The current state is not one bug but three different situations:

### 1. `npglpreg`

Status:

1. Already has explicit outer restart handling.
2. Already stores restart diagnostics.
3. Missing only the rich `np`-style one-line NOMAD progress stream.

Implication:

1. This is primarily an I/O ownership/detail problem, not a restart-semantics problem.

### 2. `crs()`

Status:

1. Still uses older `frscvNOMAD()` / `krscvNOMAD()` plumbing.
2. Still relies on `snomadr(..., nmulti=...)` for multistart behavior.
3. Already emits detailed old-style NOMAD output via `DISPLAY_DEGREE`.

Implication:

1. The first job is to determine whether the old path is only cosmetically wrong or semantically wrong.
2. If the user-visible information is already present but rendered badly, the lowest-risk fix is an I/O adapter, not restart-engine replacement.

### 3. `crsiv` / `crsivderiv`

Status:

1. Outer IV-stage progress exists.
2. Nested `crs()` NOMAD visibility is suppressed.

Implication:

1. This is a progress ownership problem only.
2. It should not require changing IV estimation or restart mechanics.

## Revised Goal

Restore the `np`-style NOMAD user experience in `crs` with the smallest possible internal disturbance.

That means:

1. Same one-line managed NOMAD progress line.
2. Same high-signal fields where available:
   - multistart
   - eval or iteration count
   - current degree
   - best degree so far
   - current/best objective value
3. Clear finish behavior.
4. No change to non-NOMAD output.
5. No multistart-engine rewrite unless current multistart behavior is proven wrong.

## Revised Plan

### Stage 0: Freeze Contracts Before Editing

This is mandatory.

Capture pre-change behavior for these routes:

1. `npglpreg(y ~ x, cv="degree-bandwidth", nmulti=3, ...)`
2. `crs(y ~ x, cv="nomad", nmulti=3, ...)`
3. `crsiv(..., display.nomad.progress=TRUE, nmulti=...)`
4. `crsivderiv(..., display.nomad.progress=TRUE, nmulti=...)`

For each route, record:

1. live progress transcript via the existing shadow renderer
2. final fitted object summary fields
3. restart-related returned-object fields, if any
4. numerical outputs needed for parity

Add focused tests that encode current solver semantics separately from desired I/O semantics.

Rule:

1. No behavioral edits before these fixtures exist.
2. Save the raw captured traces under a dated `/tmp` root before the first code change.
3. Include at least one fixed-seed objective/result parity artifact for each touched route.

### Stage 0.5: Predeclare Hard Stop Criteria

Immediate rollback for the active tranche if any of the following occur:

1. fixed-seed selected degree, segments, bandwidths, or chosen restart drift for the touched route,
2. fixed-seed objective value drifts beyond floating-point noise,
3. non-NOMAD output changes for the touched public entry point,
4. nested progress lines duplicate, disappear, or fail to clear,
5. a new shared helper changes untouched routes.

### Stage 1: `npglpreg` I/O Only

Scope:

1. No restart logic changes.
2. No start generation changes.
3. No objective-function changes.
4. No `snomadr` changes.
5. No shared-helper extraction before route-local proof.

Implementation target:

1. First implement the formatter locally inside `R/np.regression.glp.R`.
2. Feed it from `glpcvNOMAD()` using data it already has:
   - outer multistart index
   - objective evaluations seen in the callback
   - current degree candidate
   - best-so-far degree/objective
3. Only after this passes validation, decide whether any tiny formatter belongs in `R/progress.R`.

Design constraint:

1. `glpcvNOMAD()` should remain the owner of the visible NOMAD line.
2. The existing coarse `Calling NOMAD` / `fv = ...` status updates should be replaced, not layered on top.

Exit criteria:

1. `npglpreg` shows one managed NOMAD line with `multistart`, `eval`, `deg`, `best`, and `fv`.
2. Returned objects are numerically unchanged.
3. Existing restart bookkeeping is unchanged.

### Stage 2: `crs()` I/O Adapter Only

This is the most important correction to the first plan.

Do not rewrite `frscvNOMAD()` / `krscvNOMAD()` multistart yet.

Instead:

1. Preserve current `snomadr(..., nmulti=...)` execution semantics.
2. Intercept the old NOMAD callback/detail stream and re-render it through the same managed-line formatter used for `npglpreg`.
3. Normalize the old `DISPLAY_DEGREE` information into the new one-line renderer.

Why this is safer:

1. The old path already appears to know current degree, segments, and CV value.
2. If we can re-own that information at the R progress layer, we may get the desired UX without touching solver semantics at all.
3. That sharply reduces regression risk.

Only if this adapter proves impossible should we reopen restart-engine replacement.

Additional constraint:

1. Treat `frscvNOMAD()` and `krscvNOMAD()` as separate tranches if their detail streams differ materially.
2. Do not “fix both because they look similar” without first proving that they expose the same recoverable fields.

Exit criteria:

1. `crs(..., cv="nomad")` shows one managed NOMAD line.
2. Underlying selected parameters remain unchanged on fixed seeds.
3. No restart/mechanics drift is introduced.

### Stage 3: `crsiv` / `crsivderiv` Ownership Handoff Only

Scope:

1. Keep outer IV-stage progress.
2. Stop muting nested NOMAD visibility when nested `crs()` calls are in a NOMAD route.
3. Resume outer stage line after nested fit completion.

Important constraint:

1. Do not change IV iteration logic.
2. Do not change the `nmulti` routing rules in `.crsiv_prepare_dot_args()` during this tranche.
3. Do not redesign the outer IV-stage wording in the same tranche.

Exit criteria:

1. IV routes remain readable.
2. Nested NOMAD detail becomes visible.
3. No numerical drift.

### Stage 4: Multistart Semantics Audit, Not Rewrite

Only after Stages 1-3 are complete and stable:

1. audit whether `crs()` old `snomadr(..., nmulti=...)` semantics materially differ from the explicit outer-restart model in ways that matter
2. compare:
   - chosen minima
   - restart starts
   - number of evaluations
   - reproducibility under fixed seeds
3. decide whether rewrite is justified

Default outcome:

1. if behavior is acceptable, do not rewrite it near release
2. if behavior is not acceptable, isolate that as a separate post-I/O tranche

## Revised Implementation Strategy

### Shared Helper Strategy

The shared code, if any, should be limited to formatting and ownership, not restart mechanics.

Good candidates after route-local proof:

1. a progress-detail formatter in `R/progress.R`
2. small route-specific adapters in:
   - `R/np.regression.glp.R`
   - `R/frscvNOMAD.R`
   - `R/krscvNOMAD.R`

Bad candidate for the current release window:

1. centralizing all multistart generation into one new architecture before proving necessity
2. extracting route-general helpers for elegance rather than risk reduction

### Test Strategy

Add or strengthen tests in this order:

1. `npglpreg` progress contract tests
2. `crs()` NOMAD progress contract tests
3. `crsiv` / `crsivderiv` progress handoff tests
4. fixed-seed numerical parity tests for each touched route
5. explicit non-NOMAD unchanged-output tests for each touched entry point

Tests must separately assert:

1. visible output contract
2. restart bookkeeping contract
3. numerical parity contract
4. finish-clear contract

That separation is important because it lets us improve I/O without accidentally relaxing solver expectations.

## Revised Risk Register

### Low-risk changes

1. replacing coarse status text with richer managed-line text in `npglpreg`
2. adding progress-formatting helpers
3. handoff logic between outer IV-stage lines and nested NOMAD lines

### Medium-risk changes

1. re-owning `frscvNOMAD()` / `krscvNOMAD()` old NOMAD detail stream at the R layer

### High-risk changes

1. replacing `frscvNOMAD()` / `krscvNOMAD()` `nmulti` execution model
2. changing `snomadr()` multistart semantics
3. changing IV loop `nmulti` routing while also changing I/O

Highest-standard recommendation:

1. Do not take the high-risk path unless the lower-risk adapter path is proven insufficient.

## Commit Strategy

Use ultra-small checkpoints.

Recommended granularity:

1. contract tests only for `npglpreg`
2. `npglpreg` route-local I/O implementation
3. contract tests only for the first `crs()` NOMAD helper path touched
4. first `crs()` NOMAD helper I/O implementation
5. second `crs()` NOMAD helper only if needed
6. `crsiv` handoff only
7. `crsivderiv` handoff only

This is intentionally slower and safer than a combined patch.

## What We Should Not Optimize For

In this release-candidate window, do not optimize for:

1. minimal code duplication,
2. elegant unification of all NOMAD routes,
3. deep helper refactoring,
4. architecture cleanliness at the expense of proof.

We should optimize for:

1. exact behavioral restoration,
2. surgical diffs,
3. easy rollback,
4. route-by-route proof.

## Revised Success Definition

The work is successful if:

1. `npglpreg`, `crs()`, `crsiv`, and `crsivderiv` all show the intended single-line NOMAD visibility where NOMAD is active
2. non-NOMAD output is unchanged
3. restart semantics are unchanged unless separately audited and intentionally changed
4. fixed-seed numerical outputs do not drift
5. no unrelated doc/worktree churn is pulled into the implementation commits

## Practical Recommendation

If implementation starts now, the safest order is:

1. freeze behavior and add progress/parity fixtures
2. fix `npglpreg` I/O only
3. fix one `crs()` NOMAD helper path by adapter, not restart rewrite
4. fix the second `crs()` NOMAD helper path only if needed
5. fix `crsiv` handoff only
6. fix `crsivderiv` handoff only
7. postpone any multistart-engine unification until after release unless a clear semantic bug is demonstrated

This revised plan has the highest probability of restoring the desired behavior with the lowest probability of regression or collateral damage.
