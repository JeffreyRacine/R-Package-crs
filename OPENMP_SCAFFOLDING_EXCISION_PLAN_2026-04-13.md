# OpenMP Scaffolding Excision Plan

Date: `2026-04-13`

Repo: `/Users/jracine/Development/crs`

## Engineering Judgment

I am not 100% comfortable with the first draft of this plan.

From a highest-engineering-standards perspective, the earlier version was still
too willing to touch vendored NOMAD surfaces that do not need to move in order
to achieve the goal safely.

The most important criticism is this:

1. the defect is a package-contract leak,
2. so the first fix must seal the package contract at the narrowest possible
   interception points,
3. and only after that should we consider removing residual package-owned
   scaffolding or cleaning up user-facing references.

Anything broader than that increases regression risk without increasing the odds
of success proportionally.

## What Was Too Risky In The First Draft

### 1. It was still too broad in vendored surfaces

The prior draft named multiple vendored NOMAD help/attribute files as likely
edit targets. That is not the safest first move.

Why this is risky:

1. those files are upstream/vendor text surfaces,
2. changing them increases fork delta,
3. they are not required to stop the broken load path,
4. they are not required to make shipped `crs` truthful at runtime.

Highest-standard correction:

- do not edit vendored help/attribute text in the first code tranche.

### 2. It did not freeze the full shipped OpenMP contract tightly enough

The prior draft focused on `COOP_MADS_OPTIMIZATION` and
`PSD_MADS_OPTIMIZATION`, but it did not elevate `NB_THREADS_PARALLEL_EVAL > 1`
to first-class excision scope.

That is too loose.

The archived investigation already established that shipped `crs` should not be
treated as supporting OpenMP parallel evaluation beyond the safe default.

Highest-standard correction:

- treat the entire experimental OpenMP surface as unsupported in shipped `crs`:
  - `COOP_MADS_OPTIMIZATION`
  - `PSD_MADS_OPTIMIZATION`
  - `NB_THREADS_PARALLEL_EVAL > 1`

### 3. It proposed removing the `match` workaround too early

That is not a safe first-tranche edit.

Even after the experimental algorithms are sealed off, ambient `_OPENMP`
environments may still cause vendored NOMAD files to include `<omp.h>` through
upstream conditionals unrelated to `COOPMads` / `PSDMads`.

Highest-standard correction:

- do not remove the `match` workaround in `src/snomadr.h` in the first
  checkpoint.
- treat it as a compatibility shim until proven unnecessary.

### 4. It relied too much on existing sentinels

The earlier draft named good installed sentinels, but it did not require a new
targeted regression test for the exact contract we are changing.

That leaves a gap.

Highest-standard correction:

- add a dedicated negative-contract test file for unsupported OpenMP options.

### 5. It was too generous about claim strength if ambient `_OPENMP` cannot be reproduced locally

Structural proof is valuable, but if we cannot locally recreate the exact
ambient-`_OPENMP` environment that produced the failure, we must not overclaim.

Highest-standard correction:

- claim tiers must distinguish:
  - `structurally sealed and ordinary builds validated`
  - from
  - `ambient-OpenMP repro reclosed empirically`

## Revised Objective

Achieve the goal with the smallest credible code delta:

1. explicitly freeze shipped `crs` as not supporting the experimental OpenMP
   routes,
2. seal runtime admission and dispatch against ambient `_OPENMP`,
3. remove the package-owned experimental build hook we added,
4. add targeted tests that lock the contract in place,
5. avoid broad vendored churn.

## Provenance

The relevant history remains:

1. `e872726` `Upgrade crs NOMAD bridge to embedded NOMAD 4.5.0`
   - introduced the vendored NOMAD 4.5 `_OPENMP` surface.
2. `b0871f2` `NOMAD: extend tuning frontier, add OpenMP findings, and harden fr defaults`
   - added the repo-owned experimental build hook:
     `CRS_EXPERIMENTAL_OPENMP=1`
   - added the `match` workaround in `src/snomadr.h`.
3. `3179d99` `Use portable static Makevars manifest`
   - preserved exclusion of `COOPMads` / `PSDMads`,
   - but no longer preserved a clean package-owned boundary from ambient
     `_OPENMP`-gated dispatch.

## Revised Strategy

This plan follows a strict minimal-diff rule:

1. seal contract first,
2. remove dead scaffolding second,
3. clean package-owned references third,
4. do not edit broad vendored surfaces unless a hard gate proves it is required.

## Tranche Structure

### Tranche 0: Baseline And Artifact Lock

Artifact root:

- `/Users/jracine/Development/tmp/crs_openmp_scaffolding_excision_20260413`

Save:

1. `git status`
2. targeted inventories for:
   - `CRS_EXPERIMENTAL_OPENMP`
   - `_OPENMP`
   - `COOP_MADS_OPTIMIZATION`
   - `PSD_MADS_OPTIMIZATION`
   - `NB_THREADS_PARALLEL_EVAL`
3. baseline source-install log
4. baseline `library(crs)` load log
5. baseline installed `nm -u` / `otool -L` evidence for `crs.so`

No code edits yet.

### Tranche 1: Freeze The Shipped OpenMP Contract Explicitly

Objective:

- encode the package contract directly in the build, instead of relying on
  ambient `_OPENMP` behavior.

Preferred design:

1. introduce one package-owned compile-time contract macro in the generated
   build flags, for example:
   - `CRS_NOMAD_EXPERIMENTAL_OPENMP=0`
2. generate it from:
   - `tools/nomad/generate_makevars_manifest.R`
   - `src/Makevars`
3. use that macro as the single source of truth for whether experimental
   OpenMP-only algorithm routes are available in shipped `crs`.

Important rule:

- do not use raw `_OPENMP` alone as the shipped package contract anymore.

Files in scope:

1. `tools/nomad/generate_makevars_manifest.R`
2. `src/Makevars`

This tranche should not yet remove the old compatibility shim in
`src/snomadr.h`.

### Tranche 2: Seal Admission And Dispatch At The Narrowest Runtime Points

Objective:

- ensure unsupported OpenMP surfaces are unreachable even if the build
  environment defines `_OPENMP`.

Only touch the two narrowest vendored interception points:

1. `src/nomad4_src/src/Algos/MainStep.cpp`
2. `src/nomad4_src/src/Param/RunParameters.cpp`

Required behavior after this tranche:

1. `COOP_MADS_OPTIMIZATION` is unavailable in shipped `crs`
2. `PSD_MADS_OPTIMIZATION` is unavailable in shipped `crs`
3. `NB_THREADS_PARALLEL_EVAL > 1` is unavailable in shipped `crs`
4. ordinary builds continue working
5. ambient `_OPENMP` cannot reopen the missing-symbol path

Design constraints:

1. guard `COOPMads` / `PSDMads` includes and dispatch on the package-owned
   contract macro, not raw `_OPENMP`
2. fail fast with package-truthful errors for unsupported options
3. do not edit `COOPMads` / `PSDMads` implementation files themselves
4. do not delete vendored directories
5. do not touch broad vendored lock/evaluator/cache/OpenMP code here

### Tranche 3: Add Targeted Regression Tests Before Cleanup

Objective:

- lock the new contract in place with tests that directly cover the defect and
  the intended shipped behavior.

Add a new focused test file, for example:

- `tests/testthat/test-snomadr-openmp-contract.R`

Required test coverage:

1. `snomadr(..., opts = list("COOP_MADS_OPTIMIZATION" = "true"))`
   fails with a shipped-`crs` unsupported message
2. `snomadr(..., opts = list("PSD_MADS_OPTIMIZATION" = "true"))`
   fails with a shipped-`crs` unsupported message
3. `snomadr(..., opts = list("NB_THREADS_PARALLEL_EVAL" = "2"))`
   fails with a shipped-`crs` unsupported message
4. a standard supported `snomadr` call still works

Why this tranche comes before broader cleanup:

- once these tests exist, any later cleanup that reopens the leak is much easier
  to catch immediately.

### Tranche 4: Remove The Repo-Owned Experimental Build Hook

Only after Tranches 1-3 are green:

1. remove the dead `CRS_EXPERIMENTAL_OPENMP` hook from the generator/build
   contract
2. remove stale comments that advertise a local experimental OpenMP install path
3. keep the `PSDMadsMegaIteration_prev.cpp` exclusion untouched

This tranche is now low-risk because the real contract is already enforced by
the package-owned macro and the runtime interception points.

### Tranche 5: Clean Only Package-Owned User-Facing References

After the code contract is sealed and tested:

1. update package-owned documentation/reference surfaces that would otherwise
   mislead users about OpenMP availability
2. prefer a narrow package-owned clarification over broad edits to vendored
   NOMAD text

Preferred surfaces:

1. `inst/nomad/NOMAD_4_5_0_OPTIONS_REFERENCE.md`
2. package-owned comments or notes in build scripts if needed

Do not edit vendored NOMAD attribute-definition text in the first cleanup pass
unless a concrete user-facing path still exposes it.

## What We Will Not Do In This Campaign

These are explicit non-goals:

1. no deletion of vendored `COOPMads` / `PSDMads` directories
2. no broad `_OPENMP` purge across vendored NOMAD
3. no edits to upstream NOMAD CMake files
4. no removal of the `match` workaround in `src/snomadr.h` unless later proof
   shows it is dead and harmless to remove
5. no edits to user-dirty release files:
   - `CHANGELOG`
   - `DESCRIPTION`
   - `crsver`

## Validation Gates

No positive checkpoint claim without all gates for the active tranche.

### Gate A: Baseline Installed Build

1. `R CMD INSTALL` from source into a scratch library
2. `library(crs)` from that installed library

### Gate B: Static Manifest/Contract Audit

Confirm after regeneration:

1. `CRS_EXPERIMENTAL_OPENMP` is gone from live build files
2. the package-owned contract macro is present and explicit
3. `COOPMads` / `PSDMads` remain absent from `NOMAD4_OBJECTS`
4. `PSDMadsMegaIteration_prev.cpp` remains excluded

### Gate C: New Targeted Regression Tests

Run the new OpenMP-contract test file added in Tranche 3.

### Gate D: Existing Simple Installed Sentinel

Run:

- `tests/testthat/test-snomadr-args.R`

### Gate E: Representative Installed NOMAD Sentinel

Run:

- `demo/nomad_kernel.R`

If needed for speed, use a smaller installed harness first, but do not close the
campaign without one representative installed NOMAD route.

### Gate F: Binary Surface Audit

On the installed `crs.so`:

1. inspect unresolved symbols with `nm -u`
2. confirm no unexpected `PSDMads` / `COOPMads` unresolved references appear

### Gate G: Ambient-`_OPENMP` Claim Discipline

If the host can reproduce an ambient-`_OPENMP` source build safely:

1. run it
2. require install/load success before claiming empirical closure

If the host cannot reproduce it cleanly:

1. save structural proof that admission and dispatch no longer depend on raw
   `_OPENMP`
2. report the result as:
   - `structurally sealed, ordinary builds validated`
   not:
   - `ambient-OpenMP repro fully reclosed`

## File Scope For The First Real Code Checkpoint

The first checkpoint should stay as close as possible to this file set:

1. `tools/nomad/generate_makevars_manifest.R`
2. `src/Makevars`
3. `src/nomad4_src/src/Algos/MainStep.cpp`
4. `src/nomad4_src/src/Param/RunParameters.cpp`
5. one new focused test file under `tests/testthat`

Anything beyond that should be treated as escalation and justified before keep.

## Keep / Drop Rule

Keep the tranche only if all are true:

1. ordinary source install/load remains green
2. unsupported OpenMP surfaces now fail fast and truthfully
3. the exact missing-symbol leak cannot reopen through ambient `_OPENMP`
   structurally
4. standard supported NOMAD routes still pass installed sentinels
5. the code delta remains narrow and reviewable

Drop or revert immediately if any of the following happen:

1. the patch expands into broad vendored OpenMP cleanup
2. ordinary installed behavior regresses
3. we need to touch user-dirty release files to complete the fix
4. the new negative-contract tests cannot be made reliable

## Recommended Execution Order

To maximize odds of success and minimize collateral damage, execute in this
exact order:

1. add the package-owned contract macro to the generated build contract
2. seal `MainStep.cpp` and `RunParameters.cpp`
3. add the targeted regression tests
4. run baseline install/load plus tests plus one representative NOMAD route
5. only then remove the dead `CRS_EXPERIMENTAL_OPENMP` hook/comments
6. only then clean narrow package-owned user-facing references

This order is the safest because it closes the defect first, proves the new
contract second, and defers all cosmetic or documentation cleanup until after
the runtime behavior is re-earned.
