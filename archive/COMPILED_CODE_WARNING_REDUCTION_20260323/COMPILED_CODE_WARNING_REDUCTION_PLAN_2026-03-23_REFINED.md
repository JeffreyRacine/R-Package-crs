# Compiled-Code Warning Reduction Plan (Refined For RC Safety)

Date: 2026-03-23
Repo: `/Users/jracine/Development/crs`
Priority: release-candidate safety
Scope: compiled-code warning reduction only

## Governing Principle

NOMAD is the heart of the package.

Therefore this tranche must optimize for:

1. maximum solver safety,
2. minimum source disturbance,
3. measurable warning reduction only when the reduction is proven harmless,
4. immediate stop when the remaining work ceases to be obviously low-risk.

This is not a cleanup campaign and not a vendor-modernization campaign.

## Stronger Critique Of The Prior Plan

The revised plan was better, but it can still be improved in four important ways.

### 1. It still framed success too much in terms of warning reduction

For RC work, the real success criterion is:

- preserve exact package behavior,
- while safely removing only warning-causing code that is proven irrelevant to package runtime.

That sounds similar, but it leads to stricter decisions:

- if the next warning source is runtime-essential or even plausibly runtime-essential, stop;
- do not keep pushing just because some warning remains.

### 2. It did not lock solver identity tightly enough

For NOMAD, "no behavioral drift" should be made more concrete.

We should freeze, per checkpoint:

- same optimization engine selection,
- same parameter acceptance behavior,
- same stop reason classes,
- same selected model structure where deterministic,
- same objective values to within explained floating-point noise only.

This is stricter than ordinary smoke success.

### 3. It did not distinguish build-graph proof from runtime proof sharply enough

We need two separate proofs before any exclusion keep:

1. build-graph proof:
   - removing the translation unit does not break compilation or linkage.
2. runtime proof:
   - package-owned routes do not depend on the removed unit in execution.

Both must hold. A successful build alone is not enough.

### 4. It did not isolate OpenMP strongly enough

The current OpenMP scaffolding should be treated as a separate issue.

Important finding from the current tree:

- the experimental OpenMP path is only enabled when `CRS_EXPERIMENTAL_OPENMP=1`,
- the compiled-code warning objects reported by `R CMD check` are not the OpenMP-gated `PSDMads`/`COOPMads` objects,
- therefore the compiled-code warning is not being triggered by the dead OpenMP scaffolding.

That means:

- do not mix OpenMP cleanup with this tranche,
- do not touch `CRS_EXPERIMENTAL_OPENMP` logic in the first checkpoint,
- do not attribute compiled-code warning reduction to OpenMP cleanup.

At most, we may note that the OpenMP scaffolding contributes to the separate GNU Makevars warning surface, but not to the compiled-code warning under discussion.

## Explicit OpenMP Finding

Based on the current sources:

- [`/Users/jracine/Development/crs/src/Makevars`]( /Users/jracine/Development/crs/src/Makevars) only activates OpenMP flags when `CRS_EXPERIMENTAL_OPENMP=1`
- otherwise it excludes `PSDMads` and `COOPMads` sources from the build
- the compiled-code warning object list you posted is dominated by:
  - core NOMAD runtime objects,
  - SGTELIB runtime objects,
  - `NomadStdCInterface.o`
- it is not dominated by `PSDMads` or `COOPMads`

Engineering conclusion:

- the dead OpenMP scaffolding did not trigger the compiled-code warning we are trying to reduce
- it should be handled, if at all, in a separate tranche

## Refined North Star

The safest path is:

1. prove which warning-causing translation units are genuinely non-runtime baggage,
2. exclude only those units, one at a time,
3. revalidate solver identity after each exclusion,
4. stop if the remainder is in runtime-essential code.

That gives the highest probability of useful improvement with the lowest probability of collateral damage.

## Hard Safety Rules

### Rule 1: No vendor source edits until exclusion options are exhausted

The first implementation checkpoint may edit only:

- [`/Users/jracine/Development/crs/src/Makevars`]( /Users/jracine/Development/crs/src/Makevars)

No vendored `.cpp` or `.hpp` files may be edited before:

1. evidence pack is complete,
2. exclusion-only checkpoint has been attempted,
3. exclusion-only warning reduction has been measured,
4. we have explicitly concluded that the remaining warning is concentrated in a very small, proven output-only runtime surface.

### Rule 2: No OpenMP changes in this tranche

Do not edit:

- `CRS_EXPERIMENTAL_OPENMP`
- `PKG_CXXFLAGS += -fopenmp ...`
- `PKG_LIBS += ... -lomp`
- `PSDMads` / `COOPMads` inclusion logic

unless and until we open a separate OpenMP tranche.

### Rule 3: One translation unit per exclusion checkpoint

If excluding vendor sources, do it one file at a time.

Initial candidate set should start with the strongest case only:

- `src/nomad4_src/ext/sgtelib/src/Tests.cpp`

No bundled exclusions in one patch.

### Rule 4: Warning-family accounting is mandatory

For each checkpoint, maintain before/after accounting for:

1. `cout/cerr`
2. `printf/puts/putchar`
3. `_exit`
4. `_rand`

If a checkpoint reduces none of these in a measurable way, revert it.

### Rule 5: Solver-identity parity is a hard gate

A checkpoint fails if it changes any of the following unexpectedly on fixed-seed canonical fixtures:

- selected degree
- selected segment count
- best objective value
- stop reason shape
- warning/error behavior
- progress ownership/I-O contract

### Rule 6: Stop early is success

If the safe exclusion path yields only modest reduction but proves the remainder is runtime-essential, stopping is the correct engineering outcome.

## Refined Execution Plan

### Stage 0: Evidence Lock

Read-only.

Save under a dated `/tmp` root:

1. current `R CMD check --as-cran` log
2. current `nm -u` unresolved-symbol report
3. warning object-family inventory
4. current `Makevars`
5. current route-smoke outputs for:
   - `crs(...)`
   - `npglpreg(...)`
   - `crsiv(...)`
   - `crsivderiv(...)`

Also save:

6. a note stating that OpenMP scaffolding is currently inactive by default and is not the compiled-code warning root cause

Exit criterion:

- evidence pack complete

### Stage 1: Candidate Exclusion Proof

Read-only.

For `Tests.cpp`, prove:

1. no package-owned source includes or calls it directly
2. no remaining linked object needs symbols defined there
3. upstream includes it as library baggage rather than as runtime prerequisite for our package routes
4. excluding it should remove at least part of the `cout`/`exit` warning surface

If any of the above is ambiguous, do not exclude it.

Exit criterion:

- explicit yes/no decision on `Tests.cpp`

### Stage 2: Exclusion-Only Pilot

Edit:

- [`/Users/jracine/Development/crs/src/Makevars`]( /Users/jracine/Development/crs/src/Makevars) only

Change:

- exclude `Tests.cpp` only

Validation:

1. clean rebuild
2. tarball-first `R CMD check --as-cran`
3. unresolved-symbol diff
4. warning-family diff
5. fixed-seed canonical route parity
6. live I/O contract check on touched NOMAD routes

Keep only if:

1. no new warning/note/error
2. no solver-identity drift
3. warning-family reduction is measurable

Otherwise revert immediately.

### Stage 3: Reassess Or Stop

After Stage 2, decide explicitly:

1. Stop, if remaining warnings are runtime-essential enough that further edits are risky.
2. Continue only if remaining warning sources are:
   - few,
   - file-local,
   - output-only,
   - and proven package-executed.

This stage should produce a short decision memo, not code by default.

### Stage 4: Optional Output-Only Runtime Patch

Only if Stage 3 says continue.

Requirements:

1. one source file only
2. one symbol family only
3. output-only callsites only
4. no exception, stop, or RNG logic touched

Preferred target order:

1. isolated `std::cout`/`std::cerr` diagnostics in clearly runtime-safe files
2. only afterward, if justified, isolated `printf` paths

Do not touch `_exit` or `_rand` in this tranche by default.

## Refined Validation Standard

### Build And Check

Required for every kept checkpoint:

1. `R CMD build crs`
2. `R CMD check --as-cran crs_0.15-41.tar.gz`
3. explicit warning diff against baseline

### Solver Identity

For fixed-seed canonical runs, compare:

1. selected degree/segment
2. objective value
3. high-level stop/termination behavior
4. returned object shape where practical

### User I/O Contracts

Must remain unchanged:

1. `crs` rich one-line NOMAD progress
2. `npglpreg` rich one-line NOMAD progress
3. `crsiv` outer-IV progress ownership
4. `crsivderiv` outer-IV progress ownership

### Regression Policy

Any of the following is a hard fail:

- any numerical drift not explained by floating-point noise
- any change in route warnings/errors
- any broken or altered progress contract
- any new `--as-cran` warning/note/error

## Best-Standard Recommendation

The highest-probability, lowest-collateral-damage next move is now very clear:

1. do not touch OpenMP scaffolding in this tranche
2. complete the evidence pack
3. attempt only `Tests.cpp` exclusion in `Makevars`
4. rebuild and diff
5. stop and reassess before any vendor runtime edits

If that yields only a small reduction, that is still useful information and may be the correct stopping point near RC.
