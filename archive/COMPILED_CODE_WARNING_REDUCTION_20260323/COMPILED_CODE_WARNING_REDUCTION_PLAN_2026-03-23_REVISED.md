# Compiled-Code Warning Reduction Plan (Revised)

Date: 2026-03-23
Repo: `/Users/jracine/Development/crs`
Scope: compiled-code warning reduction only

## Why The First Draft Was Not Conservative Enough

The first draft was directionally sound, but it was still too permissive for release-candidate work in several ways.

### 1. It moved from classification to vendor patching too quickly

The earlier draft allowed a fairly direct progression from:

- classify runtime callsites,
- then patch "diagnostic-only" output.

That is still too aggressive. In vendored numerical code, many apparently diagnostic lines are entangled with:

- error-path assumptions,
- upstream debugging hooks,
- latent feature activation,
- or state/reporting relied on indirectly by surrounding code.

For RC safety, "looks diagnostic" is not enough. We need:

1. reachability proof,
2. runtime-irrelevance proof for package use,
3. one-file-at-a-time edits only after exclusion options are exhausted.

### 2. It did not define a strong enough proof burden for source exclusion

The first draft correctly identified `Tests.cpp` as a dead-weight candidate, but it did not require a strict build-graph proof before exclusion.

That is too weak. For this tranche, a file should be excluded only if we can show all of the following:

1. it is not referenced by package-owned source,
2. it is not needed by the vendored library objects that remain,
3. removing it does not alter linked object closure except by removing that file's symbols,
4. no runtime feature exercised by package smokes changes.

### 3. It did not separate "warning reduction" from "warning elimination" sharply enough

This matters because the safest near-RC outcome may be:

- reduce the warning materially,
- stop,
- and defer the risky remainder.

The earlier plan could still pull us into diminishing-return edits chasing a perfect result.

For highest engineering discipline, the plan must prefer:

- safe partial reduction with hard proof,

over:

- broader edits for a cleaner check log at disproportionate risk.

### 4. It did not explicitly require object-file and symbol-family accounting

We need more than grep-based source awareness.

For each step, we should know:

1. which object files are responsible for each warning family,
2. which warning families survive after each checkpoint,
3. whether the remaining families are in package-reachable runtime paths.

Without that accounting, we risk making source edits that feel justified but barely affect the real warning surface.

### 5. It did not make "no semantic vendor edits in early stages" explicit enough

Near RC, changing vendored runtime code should be a last resort.

The earlier draft allowed low-risk vendor silencing too early. The safer standard is:

- first prefer exclusion,
- then prefer compile-time containment,
- only then consider source edits,
- and only for one symbol family at a time.

### 6. It did not explicitly ban mixed-risk commits

For this work, the following must never be combined in one commit:

- `Makevars` exclusions,
- runtime vendor output silencing,
- error-path changes,
- RNG changes,
- C-interface exception-path changes.

Each of those is a separate risk axis.

## Revised North Star

The goal is not "make the warning disappear at any cost."

The goal is:

1. achieve the largest safely provable warning reduction,
2. with zero behavioral regression,
3. using the smallest possible edits,
4. and stop before semantics risk outweighs check-log benefit.

That is the highest-standard release-candidate objective.

## Hard Rules

### Rule 1: Read-only proof before code edits

No source edits until we have a read-only evidence pack containing:

1. current warning inventory,
2. current unresolved-symbol inventory,
3. object-file attribution by warning family,
4. candidate exclusion list with dependency reasoning,
5. package-reachable runtime route list for affected NOMAD/SGTELIB features.

### Rule 2: Exclusion before modification

If a warning source can be removed from the package build safely, prefer exclusion over source patching.

### Rule 3: No semantic vendor edits in the first checkpoint

The first implementation checkpoint may change:

- [`src/Makevars`](/Users/jracine/Development/crs/src/Makevars)

and nothing else.

No vendored `.cpp` or `.hpp` file edits are allowed in the first checkpoint.

### Rule 4: One symbol family per checkpoint after that

If later stages are needed, treat these as separate checkpoint families:

1. `std::cout` / `std::cerr`
2. `printf` / `puts` / `putchar`
3. `_exit`
4. `_rand`

Never mix these in one commit.

### Rule 5: High-risk families are opt-in, not default

`_exit` and `_rand` are deferred by default.

They should be touched only if:

1. earlier safe reductions are insufficient,
2. reachability is proven,
3. exact behavioral preservation is demonstrable,
4. we explicitly decide the residual warning justifies the risk.

## Revised Staging Plan

### Stage A: Baseline Lock

This stage is read-only.

Record:

1. `R CMD check --as-cran` full log
2. unresolved-symbol scan of `crs.so`
3. object-file listing contributing to each warning family
4. current `Makevars`
5. current runtime behavior for:
   - `crs(...)` NOMAD route
   - `npglpreg(...)` NOMAD route
   - `crsiv(...)`
   - `crsivderiv(...)`

Freeze not only outputs, but also:

- selected degree/segment,
- objective values,
- route warnings,
- progress-line ownership.

Stage exit criterion:

- evidence pack complete and archived under a dated `/tmp` root.

### Stage B: Build-Graph And Reachability Census

This stage is still read-only.

Produce a table with one row per warning-contributing object family:

- object file
- source file
- warning family (`cout/cerr`, `printf`, `exit`, `rand`)
- upstream build status
- package runtime reachability status:
  - proven reachable,
  - plausibly reachable,
  - unproven,
  - provably unreachable
- proposed action:
  - keep,
  - exclude,
  - defer,
  - investigate

This stage must also answer:

1. Are there any clearly dead translation units whose exclusion changes no runtime-linked closure?
2. Which warning families are dominated by those dead units?
3. Which remaining warning families live only in active runtime objects?

Stage exit criterion:

- table complete enough to justify or reject a first exclusion-only checkpoint.

### Stage C: Exclusion-Only Pilot

This is the first code-edit stage.

Allowed file changes:

- [`src/Makevars`](/Users/jracine/Development/crs/src/Makevars) only

Allowed action:

- exclude one pre-approved dead-weight source candidate at a time

Initial candidate set:

- `src/nomad4_src/ext/sgtelib/src/Tests.cpp`

Not allowed:

- excluding multiple candidates at once,
- any vendored source edits,
- any compile-flag changes beyond the minimum exclusion necessary.

Validation for each exclusion candidate:

1. clean rebuild
2. tarball-first `R CMD check --as-cran`
3. warning diff
4. focused runtime smokes
5. fixed-seed parity checks on canonical NOMAD routes

Keep criteria:

1. build succeeds,
2. no new warning/note/error,
3. no behavioral drift,
4. warning surface is reduced in a measurable way.

Revert criteria:

1. any drift,
2. no measurable warning reduction,
3. ambiguous linkage/runtime impact.

Stage exit decision:

- if the reduction is meaningful enough for RC, stop here;
- otherwise continue only if the residual warning is concentrated in demonstrably safe-to-silence runtime callsites.

### Stage D: Runtime Output Silencing Census

This stage is read-only again.

Only if Stage C is insufficient, produce a second classification table for the remaining warning sources, but with a stricter rubric:

For each surviving source occurrence, classify as:

1. package-unreachable
2. runtime-reachable but non-executed in package smokes
3. executed in package smokes but output-only
4. executed in package smokes and intertwined with control flow
5. unknown

No code edits in this stage.

Stage exit criterion:

- at least one surviving symbol family has a narrow, clearly output-only, package-executed source set.

If not, stop and defer.

### Stage E: Output-Only Vendor Patch Pilot

This stage is optional and must be extremely narrow.

Allowed only if Stage D proves:

1. the callsite is package-executed,
2. the callsite is output-only,
3. the callsite does not alter return values, stop conditions, or branching,
4. no build-exclusion alternative remains.

Allowed scope for one checkpoint:

- one symbol family,
- in one source file,
- with one validation pass.

Examples of acceptable shape:

- suppressing an unconditional debug `std::cout` line in one proven runtime file

Not acceptable:

- sweeping replacement across many vendor files,
- changing exception/reporting semantics globally,
- touching both C interface and SGTELIB core in the same patch.

Keep criteria:

1. warning family shrinks measurably,
2. all parity and I/O gates hold,
3. patch diff remains trivially reviewable.

Otherwise revert immediately.

### Stage F: Explicitly Deferred Families

These are not part of the default plan and require separate discussion if reached:

1. `_exit`
2. `_rand`
3. C-interface `printf` exception paths if they are part of actual runtime error handling

Default action:

- document and defer.

## Revised Validation Pack

### Required for Every Kept Checkpoint

1. clean tarball build
2. tarball-first `R CMD check --as-cran`
3. before/after warning-family diff
4. before/after unresolved-symbol diff
5. focused runtime smoke pack
6. fixed-seed parity pack

### Focused Runtime Smoke Pack

Minimum:

1. `crs(...)` NOMAD route
2. `npglpreg(...)` NOMAD route
3. `crsiv(...)`
4. `crsivderiv(...)` as far as the standing unrelated issue allows

Check not only success, but:

- selected degree/segment,
- objective values,
- warning text,
- progress ownership and formatting.

### Fixed-Seed Parity Pack

For at least one representative canonical fixture per route:

- same seed,
- same controls,
- compare returned selections and objective values before vs after.

If equality is not exact, quantify drift and treat any unexplained drift as a fail.

## Revised Commit Policy

Checkpoint commits must be:

1. single-axis,
2. independently revertible,
3. independently validated,
4. small enough to review in one sitting.

Allowed commit sequence:

1. exclusion-only commit
2. optional output-only vendor patch commit

If either sequence expands beyond that, the tranche should pause for reassessment.

## Revised Stop Conditions

Stop and discuss before further edits if any of the following occur:

1. the first exclusion-only checkpoint yields negligible warning reduction,
2. remaining warning families are concentrated in clearly runtime-essential objects,
3. any candidate patch touches control flow, exceptions, or randomness,
4. any parity drift appears,
5. the only way to continue is broad vendor surgery.

## Recommended First Checkpoint

The safest first checkpoint is now explicitly:

1. complete the Stage A/B evidence pack,
2. test a `Makevars`-only exclusion of `Tests.cpp`,
3. rebuild, recheck, and diff the warning inventory,
4. stop and reassess before any runtime vendor source patching.

That checkpoint has the highest chance of producing useful information with the lowest chance of collateral damage.
