# crs NOMAD Native API Hardening Checkpoint

Date: 2026-05-27

Scope: `crs` only, focused on the final native callable
`crs_nomad_solve()` and its public package-author contract.

This checkpoint is opened before continuing downstream `np` / `npRmpi`
native NOMAD route migration. Its purpose is to harden the foundation, not to
expand the migration surface.

## Trigger

The forensic review in:

- `issue_notes/crs_nomad_native_forensic_review_20260527.md`

identified several plausible native-API defects. Some are over-classified, but
enough are credible that the next engineering-safe move is a narrow `crs`
hardening checkpoint before more estimators are placed on top of the API.

## Original Requested Checklist

Do only this:

1. Add direct `crs_nomad_solve` tests.
2. Harden R callback protection and cleanup.
3. Reject or force `NB_THREADS_PARALLEL_EVAL=1` for R callbacks.
4. RAII-wrap NOMAD problem/result handles.
5. Fix option parsing/mirroring.
6. Revisit `Inf` output handling with explicit `OBJ` / `PB` / `EB` tests.
7. Add docs for seed, bounds-as-infinity, categorical enum behavior, and
   callback threading.

## Critical Assessment

The checklist is directionally right, but needs tightening before execution.

### 1. Tests must lead, not follow

The highest-risk mistake would be to patch the native API based on static
analysis alone and then rely on downstream `np` / `npRmpi` smoke tests as proof.
Those consumer tests are valuable, but they do not isolate the `crs` ABI
contract. The checkpoint must begin with direct `crs_nomad_solve()` tests that
can fail on the current implementation or prove an invariant before patching.

### 2. "Reject or force" is not neutral

Silently forcing `NB_THREADS_PARALLEL_EVAL=1` would violate the workspace
contract against overriding explicit user-selected optimizer options. If a
caller explicitly requests an unsafe thread count for an R callback, the correct
behavior is a precise error, not a silent rewrite. If no option is supplied,
NOMAD's default single-thread behavior is already acceptable.

### 3. R callback hardening is the only potentially delicate repair

`native_eval_r()` touches the R API inside a NOMAD callback. The hardening must
respect both R's longjmp behavior and C++ cleanup. A patch that merely adds C++
RAII objects around R API calls is insufficient if R longjmp bypasses C++
destructors. The repair needs an R-aware cleanup strategy, or the code should
avoid introducing cleanup dependencies across longjmp-prone calls.

### 4. `Inf` output handling is a semantics change unless proven

Allowing `Inf` is likely correct for constraint outputs, but this must be
tested by output type. `NaN` should remain a callback failure. `Inf` in
`PB`/`EB` constraints should be tested as an infeasibility signal, while
`OBJ = Inf` should be compared against NOMAD/snomadr behavior before changing
the contract.

### 5. ABI compatibility should be preserved

This checkpoint should not redesign the native ABI. Do not add block callbacks,
new output types, new seed fields, or public convenience fields here. Clarify
and harden the current API; leave expansion to a separate design tranche.

### 6. Shared-helper extraction is tempting but too broad

Factoring duplicated helpers between `snomadr.cpp` and `snomadr_native.cpp`
would be good long-term engineering, but it increases blast radius. This
checkpoint should only extract a helper if it is required to repair an accepted
defect safely. Otherwise, duplicate code can be cleaned later.

### 7. Downstream proof is still required after crs hardening

Even if this checkpoint edits only `crs`, the native API now underpins
validated `np` / `npRmpi` routes. After the direct `crs` tests pass, run a
compact installed downstream smoke to confirm no regression in representative
already-migrated routes. Do not edit `np` / `npRmpi` in this checkpoint.

## Reformulated Plan

## Second-Pass Refinement

After further criticism, the plan is acceptable only with the following
additional guardrails. These are intended to maximize the odds of a correct
repair while minimizing regression and collateral damage.

### Test Architecture Must Be Explicit

Direct `crs_nomad_solve()` tests need a small package-internal test bridge.
The preferred shape is:

1. add one unexported `.Call` test shim in `crs` source that exercises
   `crs_nomad_solve()` directly;
2. keep the shim narrow and clearly named as test/support code;
3. use it from `tests/testthat`;
4. avoid compiling ad hoc temporary packages during tests;
5. avoid exposing a new user-facing R API.

This is cleaner than relying only on downstream `np` / `npRmpi` routes, and
less fragile than making `testthat` dynamically build a consumer package.

### Commit/Checkpoint Granularity

If this checkpoint is implemented, commits should be route-sized and reversible:

1. direct tests and any test shim;
2. thread-option rejection for R callbacks;
3. RAII cleanup and R callback protection;
4. option parsing/mirroring;
5. `Inf` output handling, only if direct tests justify a behavior change;
6. documentation updates;
7. compact downstream smoke artifacts/control-note update.

If any commit fails its focused proof, revert that commit before proceeding.

### R Callback Hardening Must Avoid Fake Safety

Do not claim longjmp safety merely because C++ RAII objects were added.

Accepted repair evidence must show:

1. an R callback error returns a structured callback failure;
2. a second clean solve in the same R process succeeds immediately afterward;
3. no stale busy flag remains;
4. no NOMAD handles are leaked on ordinary C++ exception paths;
5. any use of R cleanup machinery is tested in the failure path it is meant to
   protect.

If a full R-longjmp-safe design requires a larger architectural change, stop
and document the required design rather than shipping partial safety language.

### Option Semantics Must Preserve Explicit User Intent

The checkpoint may reject unsafe explicit options, but must not silently rewrite
them.

Specifically:

- explicit `NB_THREADS_PARALLEL_EVAL > 1` with `CRS_NOMAD_CALLBACK_R` must fail
  with a precise message;
- explicit `NB_THREADS_PARALLEL_EVAL = 1` must remain accepted;
- absent `NB_THREADS_PARALLEL_EVAL` must remain default behavior;
- explicit caller `MAX_EVAL` must not be overwritten by compatibility mirroring;
- `MAX_BB_EVAL` mirroring remains a compatibility helper only when `MAX_EVAL`
  is absent.

### Inf Output Handling Needs A Decision Table

Before changing output validation, write down and test the intended table:

| Output value | OBJ | PB | EB |
| --- | --- | --- | --- |
| finite | accepted | accepted | accepted |
| `NaN` | callback failure | callback failure | callback failure |
| `+Inf` / `-Inf` | characterize before changing | infeasibility candidate | infeasibility candidate |

The table must be confirmed against NOMAD behavior before any acceptance claim.
If `OBJ = Inf` has ambiguous behavior, leave it unchanged or reject it with
documentation rather than broadening silently.

### Downstream Smoke Is A Regression Guard, Not Discovery

The compact downstream `np` / `npRmpi` smoke exists only to confirm that the
`crs` hardening did not disturb already migrated routes. It must not become a
new estimator migration campaign, and failures should first be interpreted as
possible `crs` regression unless a plain installed `crs` direct test contradicts
that.

### No ABI Expansion In This Checkpoint

Do not spend this checkpoint on:

- block-evaluation callbacks;
- new output enum values;
- new seed fields;
- a public R wrapper for `crs_nomad_solve()`;
- exported helper APIs;
- shared helper extraction with `snomadr.cpp`;
- legacy code retirement.

Those may be good future engineering work, but they are not part of this
stabilization checkpoint.

### Non-Goals

Do not:

- edit `np-master` or `np-npRmpi`;
- add new estimator migration routes;
- remove legacy `snomadr()` or `.np_nomad_search()` code;
- add block-evaluation callbacks;
- add new native output enum values unless a test proves the current enum
  contract is actively wrong for supported callers;
- redesign random seed semantics;
- silently rewrite explicit user options;
- run broad downstream matrices as discovery.

### Phase 1: Direct crs Native API Tests

Add focused direct coverage for `crs_nomad_solve()` in `crs`.

Required tests:

1. C callback happy path:
   - small 1D or 2D quadratic objective;
   - finite bounds;
   - deterministic seed;
   - `max_eval` convenience budget;
   - success status, finite objective, solution length, and counter sanity.

2. R callback happy path:
   - same basic objective using `CRS_NOMAD_CALLBACK_R`;
   - protected function/environment path;
   - success status and sane counters.

3. R callback error path:
   - R evaluator throws;
   - solve returns `CRS_NOMAD_CALLBACK_FAILURE`;
   - immediately running a second clean solve in the same process succeeds.
   - This is the critical "busy flag not stuck" regression test.

4. Option mirror path:
   - `MAX_BB_EVAL` without `MAX_EVAL` synthesizes compatible behavior;
   - lower/mixed-case option names are classified consistently with NOMAD's
     case-insensitive option surface if accepted by the backend.

5. Unsafe thread option:
   - `CRS_NOMAD_CALLBACK_R` plus explicit `NB_THREADS_PARALLEL_EVAL > 1`
     fails clearly before entering NOMAD evaluation;
   - `NB_THREADS_PARALLEL_EVAL = 1` succeeds.

6. `Inf` output handling:
   - `NaN` output remains a callback failure;
   - `PB` and `EB` `Inf` outputs are classified according to NOMAD semantics,
     not as generic callback corruption;
   - `OBJ = Inf` behavior is characterized and documented before any change is
     accepted.

7. Cache isolation sanity:
   - two sequential solves with different objectives in the same R process do
     not reuse stale objective values.
   - This guards a claimed forensic issue even though the embedded C interface
     already appears to reset NOMAD state.

### Phase 2: Minimal Source Repairs

Implement only repairs supported by Phase 1.

Required candidate repairs:

1. R callback protection and cleanup:
   - ensure R callback `eval_f` and environment are safe for the duration of
     the solve;
   - ensure R callback failures cannot leave `native_solver_busy` set;
   - ensure failures do not leak NOMAD handles;
   - use R-aware cleanup where R longjmp can occur.

2. Thread option guard:
   - reject explicit `NB_THREADS_PARALLEL_EVAL > 1` for
     `CRS_NOMAD_CALLBACK_R` with `CRS_NOMAD_INVALID_INPUT` and a precise
     message;
   - do not silently force or rewrite the caller's option;
   - keep `NB_THREADS_PARALLEL_EVAL = 1` accepted.

3. RAII cleanup for NOMAD handles:
   - wrap `NomadProblem` and `NomadResult` handles so all ordinary C++
     exception and early-return paths free them;
   - keep the C ABI unchanged.

4. Option parsing/mirroring:
   - trim `nomad.opt` lines before forwarding;
   - handle option-name matching case-insensitively for internal mirror logic;
   - parse integer-like `MAX_BB_EVAL` values robustly enough to match accepted
     NOMAD inputs;
   - preserve explicit option precedence.

5. `Inf` output behavior:
   - reject `NaN`;
   - stop treating every `Inf` output as a callback failure if direct tests
     show NOMAD expects `Inf` as a valid infeasibility signal for constraint
     outputs;
   - document any remaining restriction.

### Phase 3: Documentation

Update `inst/include/crs_nomad_native.h` and `man/crs_nomad_api.Rd`.

Required documentation:

1. Seed semantics:
   - `random_seed > 0` sets NOMAD `SEED`;
   - `random_seed == 0` leaves NOMAD's seed unset unless explicit `SEED` is
     provided in options;
   - generated starts use `random_seed` as currently documented;
   - callers wanting NOMAD PID/non-deterministic seeding must use explicit
     NOMAD options if supported.

2. Bounds-as-infinity:
   - use `-INFINITY` / `INFINITY` for unbounded coordinates.

3. Categorical input enum:
   - `CRS_NOMAD_INPUT_CATEGORICAL` is currently mapped to integer because
     NOMAD 4 does not provide true categorical input support;
   - package consumers that need categorical semantics must encode the
     appropriate continuous/integer search coordinate and decode in their
     callback.

4. Callback threading:
   - C callbacks must be thread-safe if future NOMAD options enable parallel
     evaluation;
   - R callbacks are single-thread only, and explicit parallel-evaluation
     requests are rejected.

5. Callback output values:
   - `NaN` is an evaluation failure;
   - accepted `Inf` behavior is stated by output type after Phase 1 proof.

### Phase 4: Validation Gates

Acceptance requires all of:

1. `crs` direct native API tests pass from an installed build.
2. Existing `snomadr()` tests still pass, including OpenMP-option rejection.
3. `R CMD check`-level package tests pass for `crs` or, if a full check is too
   costly during the checkpoint, a documented installed `testthat` run plus a
   scheduled full check before downstream migration resumes.
4. A compact downstream smoke, using installed packages only, passes for
   representative already-migrated `np` / `npRmpi` routes:
   - one serial native C-callback route;
   - one serial native R-callback / degree-search route if present;
   - one MPI native C-callback route with `nslaves=1`;
   - one MPI native C-callback route with `nslaves=2`;
   - explicit `npRmpi.quit.returned=0` where MPI is used.

### Stop Conditions

Stop for Jeffrey only if:

1. direct tests prove the current public native ABI cannot express required
   safe behavior without adding fields or changing enum semantics;
2. `Inf` output semantics require a statistical or optimizer-contract decision;
3. a hardening patch changes downstream `np` / `npRmpi` selected solutions
   outside already accepted drift rules;
4. R longjmp-safe cleanup requires a larger architectural change than this
   checkpoint can safely contain.

Otherwise, proceed through the checkpoint without expanding scope.

## Success Definition

This checkpoint is complete when `crs_nomad_solve()` has direct installed-test
coverage for C and R callback modes, unsafe R-callback threading is rejected,
cleanup is robust on ordinary failure paths, option handling matches the
documented contract, `Inf` output behavior is tested and documented, and compact
downstream smoke confirms no regression in already migrated `np` / `npRmpi`
native routes.

## Execution Closeout

Status: green / validated for this `crs`-only hardening checkpoint.

Implemented:

1. Added an unexported direct `.Call` test bridge:
   `crs_nomad_native_test_solve`.
2. Added installed direct tests for:
   - C callback happy path;
   - R callback happy path;
   - R callback error followed by a clean solve in the same process;
   - rejection of `NB_THREADS_PARALLEL_EVAL > 1` for R callbacks;
   - trimmed `nomad.opt` parsing and R-callback thread screening;
   - `NaN`, `OBJ = Inf`, and `PB`/`EB = Inf` output handling;
   - sequential native solves with distinct objective surfaces.
3. Hardened R callback execution with `R_ToplevelExec()` / silent R evaluation
   and a valid `R_MakeUnwindCont()` token for outer cleanup protection.
4. Added single-source native callback-failure marking and cleanup-safe NOMAD
   problem/result handle release.
5. Hardened option parsing:
   - case-insensitive `MAX_EVAL` / `MAX_BB_EVAL` recognition;
   - integer-like value parsing;
   - trimmed/comment-stripped `nomad.opt` lines;
   - precise rejection of unsafe R-callback parallel evaluation.
6. Documented seed behavior, infinity bounds, categorical enum behavior,
   callback threading, and output-value semantics in the installed header and
   `.Rd` page.

Important diagnosis:

- The first direct tests exposed an apparent endless `Error: bad value` loop.
  Sampling showed this came from calling `R_UnwindProtect()` with
  `R_NilValue` instead of a continuation token. The final patch uses
  `R_MakeUnwindCont()` and clears the token result cell after normal return.
- A speculative edit to NOMAD's bundled C interface was discarded before final
  proof because it was diagnostic-only and not the root cause.

Proof artifacts:

- Private install:
  `/Users/jracine/Development/tmp/crs_nomad_native_hardening_20260527/Rlib/candidate_final`
- Final install log:
  `/Users/jracine/Development/tmp/crs_nomad_native_hardening_20260527/logs/install_candidate_final.log`
- Direct native API test log:
  `/Users/jracine/Development/tmp/crs_nomad_native_hardening_20260527/logs/test_direct_native_solve_final.log`
- Legacy `snomadr()` smoke:
  `/Users/jracine/Development/tmp/crs_nomad_native_hardening_20260527/logs/snomadr_legacy_smoke_final.log`
- Compact downstream smoke:
  `/Users/jracine/Development/tmp/crs_nomad_native_hardening_20260527/logs/downstream_smoke_private_crs.log`
- MPI `nslaves=2` downstream smoke:
  `/Users/jracine/Development/tmp/crs_nomad_native_hardening_20260527/logs/downstream_smoke_nprmpi_nslaves2_private_crs.log`

Proof summary:

- Private `R CMD INSTALL` passed.
- Direct `crs_nomad_solve()` installed tests passed.
- Legacy `snomadr()` quadratic smoke passed.
- Serial `np::npudensbw()` smoke against private `crs` passed.
- `npRmpi::npudensbw()` smoke against private `crs` passed for `nslaves=1`
  and `nslaves=2`, with `npRmpi.quit()` returning `0`.
- `git diff --check` passed.

Known note:

- The direct `PB`/`EB = Inf` characterization tests print NOMAD's own
  exception diagnostics for infeasible/degenerate constraint cases, but return
  through the native status path without hanging, poisoning the next solve, or
  misclassifying the callback as a native callback failure.
