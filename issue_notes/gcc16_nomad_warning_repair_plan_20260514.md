# GCC 16 NOMAD Warning Repair Plan

Date: 2026-05-14
Package: `crs`
Target version: next CRAN repair after `0.15-43`

## Context

CRAN asked that `crs` correct compilation WARNINGs on the
`r-devel-linux-x86_64-debian-gcc` check flavor before 2026-06-04. That flavor
uses Debian testing GCC/G++ 16 prerelease compilers.

The current CRAN check log is for `crs 0.15-42`, not the local `0.15-43` tree.
The check result is otherwise clean: the only status is an installation
WARNING, and all post-install checks pass.

The significant install warnings promoted by CRAN are GCC 16
`-Warray-bounds` diagnostics in vendored NOMAD C++ code:

1. `NOMAD_4_5::SgtelibModel[0]` partly outside `unsigned char [56]`
   - warning site: `/usr/include/c++/16/bits/stl_construct.h`;
   - inlined from `std::make_shared<NOMAD::SgtelibModel>(...)`;
   - package source context:
     `src/nomad4_src/src/Algos/Mads/SgtelibSearchMethod.cpp`.
2. `NOMAD_4_5::ProgressiveBarrier[0]` partly outside `unsigned char [184]`
   - warning site: `/usr/include/c++/16/bits/stl_construct.h`;
   - inlined from `std::make_shared<NOMAD::EvalPoint>(...)`;
   - package source context:
     `src/nomad4_src/src/Eval/ProgressiveBarrier.cpp`.

The full install log also contains many pre-existing vendored NOMAD warnings,
especially `-Woverloaded-virtual` around `ArrayOfDouble::operator<` and several
signedness warnings. CRAN's summary identifies the two `-Warray-bounds`
families as the significant installation warnings to repair.

## Critique Of The First Draft

The first draft had the right non-negotiables, but it was not yet narrow enough
for the lowest-regression repair.

Problems to correct:

1. It allowed too many source-level repair shapes. Mentioning class hierarchy
   changes such as virtual-destructor edits creates unnecessary blast radius
   for a compiler-warning repair unless a real lifetime bug is proven.
2. It treated broad cross-package validation as a single large gate. That is
   necessary before release claims, but too expensive and noisy as the first
   feedback loop.
3. It did not state a strong enough "two-site first" rule. The warning evidence
   points to two `std::make_shared` construction patterns; the safest first
   candidate should touch only those patterns.
4. It did not separate "warning disappearance" from "behavior preservation."
   Both must be true, and a patch that validates numerically but does not clear
   GCC 16 is not a repair.
5. It allowed targeted warning suppression too early. Suppression can be an
   acceptable final containment for a proven false positive, but using it
   before a minimal source-equivalent candidate would reduce the chance of CRAN
   acceptance and hide future real diagnostics.
6. It did not explicitly require an abandonment rule for candidates that touch
   vendored NOMAD broadly. Broad vendored changes are the highest collateral
   risk in this package.
7. It named every native surface, but did not classify which ones are direct,
   adjacent, or release-wide. This makes the plan sound exhaustive without
   helping decide what to run first.

The reformulated plan below reduces degrees of freedom: reproduce, patch only
the two implicated construction sites if possible, prove GCC 16 cleanliness,
then expand validation in stages.

## Second Critique And Further Refinement

The reformulated plan is better, but still not something to treat as 100%
comfortable before evidence exists.

Remaining concerns:

1. Candidate A changes allocation shape from `make_shared`'s single allocation
   to `shared_ptr(new T(...))`'s object-plus-control-block allocation. That is
   semantically equivalent in ordinary C++ ownership terms, but it is not
   byte-for-byte implementation-equivalent and could matter in a hot path if
   used repeatedly.
2. The `ProgressiveBarrier` warning is emitted while compiling code inside
   `ProgressiveBarrier::init(...)`, but the diagnostic text names
   `ProgressiveBarrier[0]`, not merely `EvalPoint`. That makes it especially
   important not to assume the apparent `make_shared<EvalPoint>` line is the
   whole mechanism.
3. A source-equivalent patch that clears GCC 16 locally is still only a
   candidate until it passes baseline/candidate behavior and timing proof.
4. The safest engineering posture is to run tiny diagnostic variants in a
   detached scratch worktree first, then keep only the smallest variant that
   both clears the warning and proves behavior-neutral.

Further refinement:

- Treat all warning-localization edits as `diagnostic-only` until proven.
- Do not patch the live `crs` tree for experiments.
- Prefer a compile-matrix of tiny variants before choosing a candidate.
- Promote only one minimal source variant to the live repo after it clears GCC
  16 and has an obvious no-semantic-change rationale.

## Blast-Radius Minimization Rule

The repair must be as close as possible to a compiler-warning neutralization at
the implicated construction expression(s), not a NOMAD maintenance tranche.

Hard blast-radius rules:

1. The production source candidate may touch only the minimal implicated
   expression(s) needed to clear GCC 16. If one expression clears one warning
   family and the other expression is unrelated, do not touch the unrelated
   expression.
2. Do not edit `src/Makevars`, compiler flags, registration code,
   `src/snomadr.cpp`, R wrappers, tests, documentation, or release metadata in
   the same commit as the source candidate.
3. Do not reformat surrounding vendored NOMAD code.
4. Do not rename variables, reorder statements, adjust comments, or modernize
   style around the touched lines.
5. Do not change includes unless the selected candidate fails to compile and a
   single missing include is the proven reason.
6. Do not add helper functions, abstractions, macros, compatibility wrappers,
   or conditional compilation for the source candidate.
7. Do not broaden validation failures into opportunistic fixes. A failed
   sentinel means the candidate is rejected or re-analyzed, not that adjacent
   behavior should be patched.
8. Version/date/CHANGELOG changes are a separate post-validation metadata
   commit, not part of the compiler-warning repair commit.

Preferred commit structure if the candidate succeeds:

1. Commit 1: minimal source candidate only.
2. Commit 2: validation artifacts or retained scripts only if they are intended
   to live in the repo; otherwise keep artifacts under
   `/Users/jracine/Development/tmp`.
3. Commit 3: release metadata only after the source candidate is validated.

If a candidate requires touching more than the minimal implicated expressions,
stop and write a new risk note before proceeding.

## Absolute Safety Contract

This tranche is a native compiler-warning repair only.

Hard constraints:

1. Do not change the numerical behavior of `snomadr()`.
2. Do not change the numerical behavior of any `crs` estimator, helper, or
   native-code route.
3. Do not change the behavior of `crs::snomadr()` when called by `np` or
   `npRmpi` for any method supporting `nomad = TRUE`.
4. Preserve all supported NOMAD methods, options, defaults, multistart
   behavior, termination behavior, status codes, evaluation counts, warnings,
   errors, messages, and return-object structure.
5. No speed or timing regressions are acceptable in touched or adjacent routes.
6. Pre/post proof must demonstrate numerical equivalence of direct
   `crs::snomadr()` calls, `crs` estimator routes, and `np`/`npRmpi`
   `nomad = TRUE` routes that depend on `crs::snomadr()`.
7. Pre/post proof must cover any registered `crs` native `.C`/`.Call` surface
   touched directly or affected by build/native changes.

Interpretation:

- If a candidate changes objective values, selected parameters, fitted values,
  bandwidths, degrees, status fields, iteration/evaluation counts, return
  names, warning/error text, or timing beyond noise bounds, reject it.
- If the warning cannot be cleared without changing NOMAD semantics, stop and
  report the blocker rather than accepting a semantic workaround.

## Non-Goals

1. Do not upgrade or replace vendored NOMAD.
2. Do not disable NOMAD, SGTELIB, model search, progressive barriers, or any
   currently shipped optimizer capability.
3. Do not change `np` or `npRmpi` code in this tranche without separate
   current-turn approval.
4. Do not alter R-layer APIs, defaults, examples, vignettes, or documentation
   semantics except minimal release metadata after the native repair validates.
5. Do not clean unrelated NOMAD warnings unless they become CRAN-blocking in
   the same GCC 16 proof.

## Forbidden Repair Shapes

Reject any candidate that does any of the following:

1. Removes NOMAD source files or object files from the build to silence GCC.
2. Changes `snomadr()` option translation, bound handling, constraints,
   callback behavior, seeds, randomization, stopping criteria, cleanup order,
   or result extraction.
3. Changes `nomad = TRUE` routing in `crs`, `np`, or `npRmpi`.
4. Converts a deterministic failure into a catch-and-continue path.
5. Adds broad warning suppression such as global `-Wno-array-bounds`.
6. Suppresses warnings for package-owned interface code such as
   `src/snomadr.cpp` or `src/crs_init.c`.
7. Edits NOMAD class hierarchy, object ownership, or destructor contracts
   unless evidence proves a real lifetime defect and Jeffrey explicitly accepts
   that larger risk axis.
8. Uses local `R CMD check` success as proof of resolution. The issue is
   compiler-flavor-specific and requires GCC 16 or refreshed CRAN-flavor
   evidence.

## Preferred Repair

The preferred first production candidate, if evidence supports it, is a
two-site source-equivalent construction rewrite. It is not automatically
accepted merely because it is small.

Hypothesis:

- GCC 16 is warning on libstdc++'s in-place `std::make_shared` control-block
  allocation/destruction analysis for these NOMAD objects.
- Replacing the warning-triggering `std::make_shared<T>(...)` calls with
  immediate `std::shared_ptr<T>(new T(...))` may preserve semantics while
  avoiding the in-place allocation pattern that triggers the warning.

Candidate A may touch only:

1. `src/nomad4_src/src/Algos/Mads/SgtelibSearchMethod.cpp`
   - replace the `_modelAlgo = std::make_shared<NOMAD::SgtelibModel>(...)`
     construction with equivalent `std::shared_ptr<NOMAD::SgtelibModel>(
     new NOMAD::SgtelibModel(...))`.
2. `src/nomad4_src/src/Eval/ProgressiveBarrier.cpp`
   - replace only the two warning-triggering `std::make_shared<NOMAD::EvalPoint>(
     evalPoint.makeSubSpacePointFromFixed(fixedVariable))` constructions with
     equivalent explicit construction.

Candidate A must not:

- change surrounding conditions, loops, cache calls, push order, or comments
  except a short note if needed;
- introduce helper functions;
- alter ownership beyond the allocation form used to create the same
  `std::shared_ptr` target object;
- touch `src/snomadr.cpp`, `src/Makevars`, or package R code.
- combine both file edits unless the diagnostic matrix shows both are required
  to clear the corresponding warning families.

If Candidate A clears GCC 16 warnings and passes staged equivalence, stop
there. Do not continue polishing vendored NOMAD.

## Diagnostic Variant Matrix

Before promoting Candidate A, test variants in a detached diagnostic worktree
under the artifact root. Label every variant `diagnostic-only`.

Required variants:

1. Baseline current `0.15-43` tarball, no patch.
2. Variant A1: replace only `_modelAlgo` construction in
   `SgtelibSearchMethod.cpp`.
3. Variant A2: replace only the two `EvalPoint` constructions in
   `ProgressiveBarrier.cpp`.
4. Variant A3: combine A1 and A2.

Optional variants only if A1/A2/A3 do not clear the warning:

1. Variant B1: introduce an explicit local `EvalPoint` temporary before the
   existing `make_shared<EvalPoint>` call, without changing allocation form.
2. Variant B2: same as B1 plus the minimal `_modelAlgo` allocation-form change.

For each variant, record:

- exact diff;
- exact compiler command;
- whether each significant warning remains;
- whether any new significant warning appears;
- whether the variant touched only the intended lines.

Promotion rule:

- promote the smallest variant that clears all significant GCC 16 warnings;
- if A1 or A2 alone clears only one warning family, combine only as needed;
- if one warning family is already absent in local `0.15-43`, do not patch its
  associated expression merely because CRAN's `0.15-42` log showed it;
- if no source-equivalent variant clears the warning, abandon source changes
  and move to last-resort containment analysis.

## Last-Resort Containment

Targeted compiler diagnostic containment is allowed only if Candidate A fails
to clear the warning or a source-level equivalent repair is proven impossible.

Rules:

1. Suppression must be limited to vendored NOMAD translation units and the
   specific GCC warning family.
2. It must not apply to package-owned interface code.
3. It must be documented as containment for a GCC 16 false positive, with the
   evidence that object lifetimes and runtime behavior are unchanged.
4. It must pass the same numerical and timing gates as a source repair.

This is a fallback, not the preferred path.

## Artifact Root

Use persistent scratch rooted under:

```text
/Users/jracine/Development/tmp/crs_gcc16_nomad_warning_20260514
```

Minimum retained layout:

```text
baseline/
candidate_A/
gcc16_logs/
parity/
timing/
release_gate/
SUMMARY.md
```

Do not use system `/tmp` for retained proof.

## Stage 0: Reproduce And Localize

Before editing:

1. Confirm clean tree:

   ```bash
   cd /Users/jracine/Development/crs
   git status --short --branch
   git rev-parse HEAD
   ```

2. Build the current local `0.15-43` tree from a clean tarball.
3. Reproduce or attempt to reproduce the GCC 16 `-Warray-bounds` warnings.
4. Save the exact compiler commands for:
   - `SgtelibSearchMethod.cpp`;
   - `ProgressiveBarrier.cpp`.
5. Confirm whether the warnings remain in local `0.15-43`.
6. Write a short localization note in the artifact `SUMMARY.md` before
   patching.

If GCC 16 is unavailable locally, record that precisely and use the closest
container or CRAN-flavor evidence available. Do not claim local reproduction if
the environment is only Apple clang or non-GCC-16.

Clean-worktree rule:

- Because this plan file itself may be untracked while being drafted, do not
  begin source experiments in the live repo until the planning note is either
  committed, stashed, or copied into the artifact root.
- Exploratory variants must be made in a detached worktree or scratch copy
  under `/Users/jracine/Development/tmp/crs_gcc16_nomad_warning_20260514`.
- The live repo receives only the selected candidate patch after diagnostic
  evidence exists.

## Stage 1: Candidate A Compile Proof

After the selected minimal variant is promoted from diagnostic-only to a
candidate patch:

1. Rebuild from a clean tarball.
2. Confirm that the two significant `-Warray-bounds` warnings are absent under
   GCC 16 or equivalent CRAN-flavor evidence.
3. Confirm no new significant compiler warning appears.
4. Confirm `R CMD check` install status no longer reports the significant
   warnings.
5. Confirm the source diff is limited to the promoted variant plus any minimal
   source-only change. Release metadata must not be included in this candidate
   diff.

If the warning remains, abandon the candidate. Do not broaden the patch without
writing a new risk note and returning to the diagnostic variant matrix.

## Stage 2: Direct `crs` Parity

Use two private libraries:

- baseline `crs` built before the patch;
- candidate `crs` built after the patch.

Run identical scripts with fixed seeds and identical inputs. Compare serialized
outputs mechanically.

Required direct `snomadr()` sentinels:

1. Single-objective smooth unconstrained problem.
2. Bound-constrained problem.
3. Nonlinear constrained problem.
4. Option-heavy path covering integer, numeric, and string option lists.
5. Multi-objective path through `smultinomadRSolve`.
6. Existing unsupported OpenMP-only option failure tests.

Acceptance:

- exact equality where the route is deterministic;
- otherwise predeclared tolerance no weaker than the existing route justifies;
- identical status/message fields, return names, dimensions, iteration counts,
  evaluation counts, and warnings/errors.

## Stage 3: `crs` Estimator And Native Surface Parity

This stage is validation only. Do not modify code in response to this stage
inside the same tranche unless the source candidate is rejected and the plan is
reopened.

Required `crs` route sentinels:

1. `crs(..., cv = "nomad")`.
2. `crsiv(..., cv = "nomad")`.
3. `crsivderiv(..., cv = "nomad")`.
4. `krscvNOMAD` and `frscvNOMAD` through ordinary user-facing calls.
5. `clsd` route that calls `snomadr()`.

Registered `.C` calls to cover:

- `RuniqueCombs`
- `gsl_bspline`
- `gsl_bspline_deriv`

Registered `.Call` calls to cover:

- `glp_model_tmm`
- `mgcv_tmm`
- `smultinomadRSolve`
- `snomadRInfo`
- `snomadRSolve`
- `crs_hat_diag`
- `crs_uniquecombs_call`
- `crs_gsl_bspline_call`
- `crs_gsl_bspline_deriv_call`

For untouched non-NOMAD native helpers, existing focused tests are acceptable if
they serialize enough output to compare baseline and candidate. For NOMAD
native calls, require direct output comparison through exported `snomadr()`
unless a direct-call test helper already exists.

## Stage 4: `np` And `npRmpi` NOMAD Parity

Validate against installed `np` and `npRmpi` without editing those repos.

Minimum `np` sentinels:

1. `npreg` / `npregbw` with `nomad = TRUE`.
2. `npcdens` / `npcdensbw` with `nomad = TRUE`.
3. `npcdist` / `npcdistbw` with `nomad = TRUE`.
4. `npqreg` with `nomad = TRUE`.
5. `npconmode` with `nomad = TRUE`.
6. `npplreg` / `npplregbw` with `nomad = TRUE`.
7. `npscoef` / `npscoefbw` with `nomad = TRUE`.
8. `npindex` / `npindexbw` with `nomad = TRUE`.
9. Existing `np` NOMAD timing-contract sentinels.

Minimum `npRmpi` sentinels:

1. Same estimator families where supported.
2. `nslaves = 1` session route for representative NOMAD calls.
3. `nslaves > 1` MPI route for representative `npreg`, `npcdens`, `npcdist`,
   `npplreg`, `npscoef`, and `npindex` NOMAD paths.
4. Existing subprocess/attach NOMAD contract tests.

Acceptance:

- selected bandwidths, degrees, objectives, fitted values, prediction surfaces,
  status fields, and evaluation counts match baseline-`crs` results;
- `npRmpi` preserves MPI execution semantics and no-serial-fallback behavior;
- no timing regression outside the mixed timing gate.

## Timing Gate

Timing must be measured after numerical parity, not used as a substitute for
it.

Rules:

1. Use repeated runs and fixed seeds.
2. Save raw per-run timings.
3. Compare median, worst-case, and absolute differences.
4. Include at least one heavier NOMAD route where allocation or optimizer
   overhead would be visible.
5. Treat persistent slowdown outside noise as candidate failure.

Expected result for Candidate A: timing indistinguishable from baseline. A
single allocation-form change in setup should not create measurable optimizer
slowdown. Because the `ProgressiveBarrier` path can allocate per cached point,
include at least one barrier/cache-heavy route so the extra-allocation risk is
measured rather than assumed away.

## Release Gate

After staged parity passes:

1. Run full `crs` source tests.
2. Build a clean tarball.
3. Install from the tarball into a private library.
4. Run the shared release gate:

   ```bash
   cd /Users/jracine/Development
   RUN_RCHK=1 ./release_protocol/run_crs_release_gate.sh
   ```

   If local `rchk` infrastructure cannot run, record the precise reason.
   `RUN_RCHK=auto` is acceptable for rehearsal only; do not claim `rchk`
   closure from an auto skip.

5. Obtain GCC 16 or refreshed CRAN-flavor evidence showing the significant
   `-Warray-bounds` warnings are gone.
6. Retain all proof artifacts under the artifact root before submission.

## Claim Vocabulary

Use precise claim tiers:

- `local signal`: reproduced, localized, or source-loaded evidence only.
- `candidate keep`: Candidate A clears the warning and passes direct `crs`
  parity plus timing gates.
- `validated`: candidate passes installed `crs`, `np`, and `npRmpi` sentinels
  plus GCC 16 warning proof.
- `release-ready`: candidate tarball passes the final release gate, required
  `rchk` proof or precisely accepted skip policy, GCC 16 warning proof, and any
  required external CRAN/win-builder checks.

Do not call the repair validated or release-ready before the corresponding
proof exists.

## Immediate Next Step

Create the artifact root, attempt a GCC 16 reproduction against the current
local `0.15-43` tarball, and write the localization note. Only then attempt
Candidate A, restricted to the two warning-triggering construction sites.
