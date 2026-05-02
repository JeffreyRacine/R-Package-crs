# rchk Native Protection Hardening Plan

Date: 2026-05-02
Package: `crs`
Target version: `0.15-43`

## Context

`crs 0.15-42` is now hosted on CRAN and regular CRAN checks are OK on the
refreshed platforms. The remaining forward-looking issue is the CRAN
"Additional issues" `rchk` report, currently shown against `crs 0.15-41`.

The report flags possible native-code protection hazards in `src/snomadr.cpp`:

1. `apply_options(...)`
   - local `SEXP names` variables remain live while calls that may allocate
     occur inside the option loops.
2. `solve_nomad4(...)`
   - local `SEXP ret` is constructed before NOMAD cleanup calls and returned
     after cleanup.

This is a native-code protection hardening issue. It is not a request to change
NOMAD behavior, option parsing, model fitting, defaults, documentation
semantics, or user-facing APIs.

## Critique Of The First Plan

The first draft had the right direction, but the implementation and validation
instructions were not yet strict enough for a CRAN-facing native-code tranche:

1. It did not state a hard clean-tree stop rule before patching.
2. It did not separate baseline and candidate installs strongly enough.
3. It said to protect `names`, but did not explicitly forbid broad rewrites of
   option parsing while touching the same code.
4. It treated performance validation too generically. A PROTECT-only patch
   should have negligible overhead; the guard should be proportionate and used
   to detect surprises, not to manufacture a performance claim.
5. It did not include GC-stress validation, which is more relevant than broad
   timing for a protection hardening change.
6. It did not explicitly classify the NSERC/SSHRC URL certificate NOTE as out
   of scope.
7. It did not state that local `R CMD check` cannot prove `rchk` resolution;
   the final `rchk` claim must come from local rchk Docker or a refreshed CRAN
   additional-issues report.

The reformulated plan below closes those gaps while keeping the blast radius as
small as possible.

## Comfort Level And Residual Risk

This plan is high confidence, but not a 100% guarantee until `rchk` itself
confirms the refreshed package. The ordinary `R CMD check` family can verify
that behavior and build surfaces remain healthy, but it cannot prove that
CRAN's static analyzer will clear every path.

Residual risks to manage explicitly:

1. `rchk` may report additional protection paths after the first sites are
   fixed.
2. A protection patch can accidentally change control flow if `UNPROTECT`
   placement is sloppy.
3. `mkcrspkg.sh` may touch generated metadata beyond the intended version/date
   fields.
4. Broad validation can create false confidence if it does not exercise the
   NOMAD option parsing paths named by `rchk`.

The constraints below are intended to maximize success odds while minimizing
regression and collateral damage.

## Goals

1. Address exactly the `rchk`-reported protection hazards.
2. Preserve behavior exactly:
   - no optimizer semantic changes,
   - no default changes,
   - no argument-handling changes,
   - no return-structure changes,
   - no warning/error/message changes except impossible crash avoidance.
3. Keep the source patch limited to the minimum native-code protection edits.
4. Validate with baseline/candidate parity, GC stress, tarball-first checks, and
   proportionate timing guards.
5. Leave unrelated CRAN notes, including acknowledgement URL certificate notes,
   out of this tranche.

## Non-Goals

1. Do not refactor `src/snomadr.cpp`.
2. Do not change NOMAD option parsing, fallback behavior, cleanup order, or
   result construction beyond the required `PROTECT`/`UNPROTECT` scopes.
3. Do not change package documentation except version/release metadata.
4. Do not combine this with URL cleanup, Rd edits, vignette edits, performance
   tuning, or release-script cleanup.
5. Do not claim the rchk issue is closed unless local rchk or CRAN's refreshed
   rchk report confirms it.

## Required Clean-Tree Rule

Before editing:

```bash
cd /Users/jracine/Development/crs
git status --short
```

If anything unrelated is dirty, stop and either commit/stash the unrelated work
or use a detached worktree rooted under `/Users/jracine/Development/tmp`.

Recommended safest workflow:

1. create a short-lived branch or detached worktree for this tranche;
2. make one source patch commit for the native protection repair;
3. make one release-metadata commit only if the source patch validates;
4. do not update CRAN staging tarballs until all validation gates pass.

## Required Version And Release Metadata Steps

1. Edit `crsver` and increment from:

   ```text
   0.15-42
   ```

   to:

   ```text
   0.15-43
   ```

2. Run:

   ```bash
   cd /Users/jracine/Development/crs
   ./mkcrspkg.sh
   ```

   Immediately inspect:

   ```bash
   git diff -- DESCRIPTION crsver CHANGELOG
   git diff --name-only
   ```

   Confirm that generated release metadata is limited and expected. If
   `mkcrspkg.sh` touches unrelated files, stop and inspect before continuing.

3. Add a new top header to `CHANGELOG`:

   ```text
   Changes from Version 0.15-42 to 0.15-43 [XX-XXX-202X]
   ```

4. Under that header, mention only the narrow hardening work:

   ```text
   * Hardened native-code protection in the NOMAD interface to address rchk
     reports, preserving existing optimizer behavior and user-facing semantics.
   ```

## Candidate Code Repair

Patch only `src/snomadr.cpp`.

Hard constraints:

1. Do not introduce helper abstractions for this patch unless explicit
   `PROTECT`/`UNPROTECT` placement becomes unreadable.
2. Prefer rchk-visible, local, explicit protection over clever C++ RAII wrappers.
3. Keep each `PROTECT` and its matching `UNPROTECT` in the same lexical block.
4. After editing, count the added `PROTECT` and `UNPROTECT` calls manually and
   verify they balance.
5. Do not put an `Rf_error`, early `return`, `break` out of the protected block,
   or new C++ exception path between `PROTECT` and its matching `UNPROTECT`.
6. If any of those constraints cannot be met cleanly, stop and redesign the
   patch before running validation.

### `apply_options(...)`

For each non-null option list (`opts_integer`, `opts_numeric`, `opts_string`):

1. Leave the existing option-loop logic unchanged.
2. Change only the `names` binding from an unprotected local to a protected
   local for the duration of that option-list loop.
3. `UNPROTECT(1)` immediately after the corresponding loop.
4. Do not add new validation, coercion, missing-name handling, fallback logic,
   or helper functions in this tranche.

Intended shape:

```c
{
  SEXP names = PROTECT(Rf_getAttrib(opts_integer, R_NamesSymbol));
  const int n = Rf_length(opts_integer);
  for (int i = 0; i < n; ++i) {
    /* existing logic unchanged */
  }
  UNPROTECT(1);
}
```

Repeat independently for `opts_numeric` and `opts_string`. Keep protection
counts local and obvious; avoid shared counters unless the final code becomes
more complex than expected.

### `solve_nomad4(...)`

Protect the return object after construction and before cleanup:

```c
SEXP ret = PROTECT(build_solution(status, message, bbe, iterations, objective, best_x));

freeNomadResult(result);
freeNomadProblem(pb);

UNPROTECT(1);
return ret;
```

Do not move cleanup, change cleanup order, change result construction, or alter
NOMAD result/problem ownership unless a reproduced failure proves the current
order is wrong.

### Explicitly Rejected Alternatives

Do not use these unless a later failed validation proves the minimal repair is
insufficient:

1. Do not rewrite option-name handling into a new vector/cache abstraction.
2. Do not change NOMAD option normalization.
3. Do not add new public or private options.
4. Do not make cleanup conditional.
5. Do not suppress or catch errors differently.
6. Do not change the return object to avoid protecting it.

## Validation Artifact Root

Retain all artifacts under:

```text
/Users/jracine/Development/tmp/crs_rchk_native_protection_20260502
```

Use subdirectories such as:

```text
baseline/
candidate/
logs/
scripts/
timing/
```

## Validation Plan

### 1. Baseline Snapshot

Before patching, record:

```bash
cd /Users/jracine/Development/crs
git status --short
git rev-parse HEAD
R CMD build .
```

Install the baseline tarball into a private library, not the user library:

```bash
R CMD INSTALL --library=/Users/jracine/Development/tmp/crs_rchk_native_protection_20260502/baseline/Rlib \
  /Users/jracine/Development/crs_0.15-42.tar.gz
```

Record `sessionInfo()` and the tarball SHA256.

### 2. Candidate Build And Install

After the version bump and protection patch:

```bash
cd /Users/jracine/Development
R CMD build crs
```

Confirm the tarball is:

```text
/Users/jracine/Development/crs_0.15-43.tar.gz
```

Install into a separate private library:

```bash
R CMD INSTALL --library=/Users/jracine/Development/tmp/crs_rchk_native_protection_20260502/candidate/Rlib \
  /Users/jracine/Development/crs_0.15-43.tar.gz
```

### 3. Focused Functional Smokes

Run the same scripts against baseline and candidate private libraries.

Minimum smokes:

1. `library(crs)` load/startup smoke.
2. One small `crs()` fit using the NOMAD-backed path.
3. One small `crscv()` or equivalent cross-validation route that exercises
   NOMAD option handling.
4. One representative example/demo already used in release preflight.

Save scripts, logs, warnings, messages, and serialized result summaries.

The NOMAD option-handling smoke must deliberately include at least one option
from each rchk-reported option list where feasible:

1. integer option,
2. numeric option,
3. string option,
4. array-style option.

If a public route cannot naturally exercise one of these categories, record
that limitation rather than adding artificial package behavior.

### 4. Parity Check

Use identical seeds, data, explicit options, and private-library paths.

Compare baseline vs candidate:

1. object classes,
2. coefficient vectors where available,
3. fitted values or predictions,
4. selected tuning parameters/options where exposed,
5. objective/CV values where exposed,
6. warnings, messages, and errors.

Acceptance:

- exact equality where deterministic objects permit it;
- documented numerical tolerance only where floating-point differences are
  expected from optimizer/platform noise;
- any structural or semantic difference is a blocker.

### 5. GC-Stress Check

Because this is a protection hardening patch, run at least one small NOMAD smoke
with forced GC stress on the candidate:

```r
gctorture(TRUE)
on.exit(gctorture(FALSE), add = TRUE)
```

Keep the dataset and iteration budget tiny. This check is for memory-protection
smoke only, not performance.

If the route is too slow under full `gctorture(TRUE)`, use the smallest
available NOMAD budget first. Only if that is still impractical, use
`gctorture2(...)` with a documented setting and state that this is a weaker
stress signal.

### 6. Performance Guardrail

Run a small paired timing guard only to detect surprising regressions.

Requirements:

1. same host,
2. same R version,
3. same BLAS/session,
4. private baseline and candidate libraries,
5. fixed seed and fixed inputs,
6. same NOMAD options and budgets.

Acceptance:

- no meaningful slowdown beyond ordinary noise;
- do not claim a speedup;
- if timings are noisy, report them as noisy and inconclusive rather than
  expanding the tranche.

### 7. Tarball-First Package Checks

Run:

```bash
cd /Users/jracine/Development
R CMD check --as-cran crs_0.15-43.tar.gz
```

Then run the package-specific release gate if available:

```bash
cd /Users/jracine/Development
./release_protocol/run_crs_release_gate.sh
```

Any new warning, error, or note is a blocker unless explicitly diagnosed and
accepted by Jeffrey.

### 8. rchk Confirmation

Preferred:

1. run local rchk Docker on `crs_0.15-43.tar.gz`, or
2. submit to win-builder/CRAN pretest and inspect the refreshed "Additional
   issues" report.

If local rchk is unavailable, the closeout must say:

- local checks cannot prove rchk clearance;
- the patch directly protects the exact reported sites;
- final rchk closure awaits a refreshed CRAN/rchk report.

Do not submit a "rchk fixed" claim to CRAN unless the refreshed rchk output is
clean or the claim is carefully phrased as "directly addresses the reported
rchk sites; awaiting refreshed rchk confirmation."

## Stop Rules

Stop and reassess before committing if any of the following happens:

1. More than `src/snomadr.cpp` plus release metadata needs editing.
2. The source diff changes more than protection placement and surrounding braces.
3. Any parity output changes.
4. The GC-stress smoke fails.
5. A package check produces a new warning/error/note.
6. The timing guard shows a surprising slowdown that cannot be explained by
   noise.
7. rchk reports a different class of issue not covered by this plan.

## Acceptance Criteria

The tranche is acceptable only if all hold:

1. Intended changed files only:
   - `crsver`,
   - `DESCRIPTION`/release metadata generated by `mkcrspkg.sh`,
   - `CHANGELOG`,
   - `src/snomadr.cpp`,
   - this plan if it is updated,
   - retained artifacts under `/Users/jracine/Development/tmp`.
2. `src/snomadr.cpp` changes are limited to protection hardening at the reported
   sites.
3. No behavior changes in focused parity checks.
4. Candidate survives the GC-stress smoke.
5. No meaningful slowdown in the paired timing guard.
6. `R CMD build` succeeds.
7. Tarball-first `R CMD check --as-cran` succeeds.
8. Release-gate artifacts are retained.
9. rchk confirmation status is honestly reported.

## Commit And Closeout Requirements

Commit message should be narrow, for example:

```text
Harden NOMAD native protection for rchk
```

Final report must state:

1. exact commit hash,
2. tarball path and SHA256,
3. files changed,
4. validation commands run,
5. baseline/candidate parity result,
6. GC-stress result,
7. performance guardrail result,
8. `R CMD check --as-cran` status,
9. rchk confirmation status,
10. any residual accepted notes or risks.
