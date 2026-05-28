# Forensic review — `crs_nomad_solve` native NOMAD API

Date: 2026-05-27
Scope: `inst/include/crs_nomad_native.h`, `src/snomadr_native.cpp`,
       registration in `src/crs_init.c`, supporting docs in
       `man/crs_nomad_api.Rd`.
Reviewer: forensic pass for latent defects, structure, engineering
quality, and computational efficiency.

## Executive summary

The native API is well-shaped at the surface: versioned header, opaque
problem/result structs with `struct_size` for forward compatibility,
reserved slots, a clear separation between C and R callback modes, and a
single registered callable. The implementation is mostly clean and
mirrors `snomadr.cpp` behaviour in the easy places.

However, the implementation has several real defects, the most serious
of which sit on the R-callback path and the exception/longjmp safety
boundary. There is also a complete absence of direct unit tests for the
new public surface — the only consumers in-tree are the
`benchmarks/nomad/*` scripts, and even those exercise `snomadr()` /
`frscvNOMAD()` / `krscvNOMAD()` rather than `crs_nomad_solve` itself.
That gap is the single highest-leverage thing to fix because nearly
every defect below would have been caught by a small suite.

Findings are tagged S1 (critical), S2 (high), S3 (medium), S4 (low /
polish).

---

## S1 — Critical defects

### S1-1. Longjmp through C++ frames leaves the global busy flag stuck

`crs_nomad_solve` guards re-entry with a process-wide
`std::atomic_flag native_solver_busy` (snomadr_native.cpp:36, 798–812).
Inside the guarded region, the R-callback path performs
`Rf_allocVector(REALSXP, nb_inputs)`, `Rf_lang2(...)`, and
`Rf_coerceVector(...)` (snomadr_native.cpp:289, 295, 316). All three
can longjmp on R-side error or allocation failure. Only the call to the
user's R function is wrapped in `R_tryEval`; the surrounding R API
calls are not.

A longjmp out of those frames will *not* execute the `try { … } catch
(…)` recovery in `crs_nomad_solve` (longjmp is not a C++ exception
under any standard-conforming implementation) and will *not* clear
`native_solver_busy`. After one such failure the process can no longer
call `crs_nomad_solve` — every subsequent call returns
`CRS_NOMAD_INVALID_INPUT` with "native NOMAD solver is already
active".

Fix: wrap the allocation/eval region in `R_UnwindProtect` (R ≥ 3.5) or
use `R_ToplevelExec` with cleanup; alternatively migrate the R-side
SEXP construction into a small helper invoked through `R_tryEval` so
that any longjmp lands inside the trap. Clear the busy flag in the
cleanup path, not after the catch.

### S1-2. R-callback SEXPs are unprotected across allocations

`native_eval_r` casts `callback->eval_f` and `callback->environment`
from `void *` to `SEXP` and uses them without `PROTECT`
(snomadr_native.cpp:275–283). The header documents that the caller
must keep them protected for the duration of the call — but the very
first thing the function does is `Rf_allocVector`, which can trigger a
GC pass. If the caller's protection is on a `PROTECT` stack that was
unwound prior to calling into NOMAD (e.g. a package author who fetched
the closure once and stored it in a long-lived `void *`), the GC will
collect it.

The fix is cheap and defensive: `R_PreserveObject` both SEXPs at the
top of `crs_nomad_solve` when `callback_mode == CRS_NOMAD_CALLBACK_R`,
and `R_ReleaseObject` on every exit. Two lines of guarded code, full
hardening of the contract.

### S1-3. R API is invoked from a thread-unsafe context with no main-thread guard

The native API is registered as a regular `R_RegisterCCallable` symbol
and the documentation does not constrain it to the R main thread. The
header explicitly says "do not call back into R, longjmp, or throw
across this C ABI" for the C callback mode — but the R callback mode
*does* call into R, on whatever thread the caller is on. R's API is
single-threaded; calling `Rf_allocVector`, `R_tryEval`,
`Rf_coerceVector` off the main thread is undefined behaviour.

This is partly mitigated because the host package that consumes the
API is typically also on the R main thread, but nothing in
`crs_nomad_solve` enforces this. NOMAD itself is multi-threaded in some
configurations; the worry isn't that the *caller* is off-thread but
that NOMAD might fan out the evaluator calls onto worker threads
internally. Verify whether `solveNomadProblem` keeps single-thread
evaluator dispatch in the current Makevars configuration, and if it
does not, the R-callback mode must serialise back to the main thread
(or be removed).

Fix: at minimum, document that `CRS_NOMAD_CALLBACK_R` is only safe
when the caller is on the R main thread. Ideally, gate the R-callback
path on `R_running_on_main_thread()` (where available) or assert the
caller is on the main thread.

### S1-4. Cross-translation-unit global state collision with `snomadr.cpp`

`snomadr.cpp` carries static globals (`thefun`, `theenv`, `routbuf`,
`rout`, `g_bbe_counter`) used by `snomadRSolve` /
`smultinomadRSolve` (snomadr.cpp:23–27). `snomadr_native.cpp` carries
its own `native_solver_busy` (line 36). The busy flag protects against
re-entrance into the *native* path only — it does *not* prevent
`snomadRSolve` and `crs_nomad_solve` from running concurrently. If
both ever fire at the same time (e.g. one R thread driving snomadr
while a worker thread reaches into the C ABI), NOMAD's internal global
state and `snomadr.cpp`'s statics will be clobbered.

Fix: either share a single process-wide lock between both back-ends, or
document the constraint and refuse simultaneous entry from the native
path while snomadr's `.Call` is in flight (track an "snomadr active"
flag in shared static storage referenced from both translation units).

### S1-5. Exception path leaks NOMAD handles

`solve_final_problem_with_starts` allocates `NomadProblem pb` and
`NomadResult nomad_result` (snomadr_native.cpp:673–688) and only
frees them on the explicit return paths. If `solveNomadProblem`,
`addNomadValParam`, `addNomadStringParam`, or any of the NOMAD C-API
helpers throws a C++ exception (these are C++ implementations under a C
veneer; nothing in the public header marks them `noexcept`), the
unwinding will bypass the `freeNomadProblem` / `freeNomadResult` calls
and leak both. The outer `catch (…)` in `crs_nomad_solve` does clear
the busy flag but does not free the handles because they're local to
the inner function.

Fix: wrap `NomadProblem` and `NomadResult` in `std::unique_ptr` with
custom deleters (one-liner with deleter functors). RAII closes every
leak path without touching the rest of the logic.

---

## S2 — High-severity defects

### S2-1. No `R_CheckUserInterrupt()` in the R-callback path

`snomadr.cpp:157` calls `R_CheckUserInterrupt()` on every evaluation
inside `r_eval_single`. `snomadr_native.cpp`'s `native_eval_r` and
`native_eval_single` do not. A long-running native solve driving an R
objective cannot be Ctrl-C'd. This is both a usability defect and an
operational one — interactive users will kill R sessions instead of
the solver.

Fix: call `R_CheckUserInterrupt()` at the top of `native_eval_r` (R
mode only; the C-callback path can't safely interrupt third-party
code). Use `R_UnwindProtect` if `R_CheckUserInterrupt` longjmps and
S1-1's fix isn't yet in place.

### S2-2. `Rf_coerceVector` inside the callback can longjmp through C++

Line 316: `rnum = PROTECT(Rf_coerceVector(result, REALSXP));`.
`coerceVector` can call into S4 method dispatch and from there into
arbitrary R code; any error/warning that reaches the top of the R
error stack will longjmp. Outside an `R_tryEval` / `R_UnwindProtect`
this is the same UB as S1-1. snomadr.cpp has the same shape but is
sheltered by being invoked from a `.Call` whose error context is set
up by R; the native API has no such shelter when called via
`R_GetCCallable`.

Fix: combine with S1-1's `R_UnwindProtect` envelope, or restrict the
acceptable return type to `REALSXP` and treat any other type as a
callback failure (faster, simpler, contract-clean).

### S2-3. `callback_evaluations` overcounts failures

`native_eval_single` increments `context->callback_evaluations` before
running the callback (snomadr_native.cpp:349), so a failed callback or
a non-finite return still increments the counter. The result struct
copies that to both `callback_evaluations` and `blackbox_evaluations`
(lines 699–700). `snomadr.cpp:198` only increments on the success
path, so `bbe` reflects accepted evaluations only.

Two callers comparing snomadr-bbe to native-bbe will see different
numbers. Fix: increment only after the callback succeeds *and* outputs
are finite, mirroring snomadr.

### S2-4. `apply_nomad_opt_file` does not trim leading whitespace

snomadr.cpp:565 trims each line before passing to `addNomadParam`.
snomadr_native.cpp:548–555 finds `first_not_of(" \t\r\n")`, uses it
only to detect comment/empty lines, then passes the **raw** line to
`addNomadParam`. NOMAD's parser is sensitive enough that a leading
space or tab will cause the option to fail, which is then converted to
`CRS_NOMAD_INVALID_INPUT` for the whole solve.

Fix: pass `line.c_str() + first` (or substr) into `addNomadParam`,
matching snomadr's behaviour.

### S2-5. `MAX_BB_EVAL` mirror is fragile

`apply_native_options` recognises `MAX_BB_EVAL` only when the value is
a clean decimal integer (snomadr_native.cpp:517–525). NOMAD itself
accepts scientific notation and surrounding whitespace via
`addNomadParam`. A caller who passes `"MAX_BB_EVAL"` with value
`"5e3"` or `" 5000 "` will not get the snomadr-compatible
auto-mirroring to `MAX_EVAL`, despite the header advertising that
behaviour ("A `MAX_BB_EVAL` option without `MAX_EVAL` follows
snomadr() compatibility by also setting `MAX_EVAL`.").

Fix: `strtol` after `strtok`-style trim of leading/trailing whitespace;
also accept `strtod` and round, mirroring `snomadr`'s
`std::lround(val)` path.

### S2-6. Option-name matching is case-sensitive

`apply_native_options` uses `std::strcmp(option->name, "MAX_EVAL")`
and `MAX_BB_EVAL` (lines 514, 517). NOMAD's internal parser is
case-insensitive. A caller passing `"max_eval"` will succeed in setting
the option but won't trigger the snomadr-compat mirror logic. Same
behavioural surprise as S2-5.

Fix: uppercase the name once, then compare. Reuse
`option_key_from_line` / `to_upper_copy` from snomadr.cpp (consider
extracting both into a shared header to avoid duplication — see
S3-2).

### S2-7. No direct test coverage for `crs_nomad_solve`

Grep across `tests/` finds zero call sites that exercise the native
callable. `test-snomadr-args.R`, `test-snomadr-cache-counters.R`, and
`test-snomadr-openmp-contract.R` cover the R-side path only. The
benchmark scripts under `benchmarks/nomad/` drive snomadr, frscvNOMAD,
krscvNOMAD — none of them go through `crs_nomad_solve`.

A public-facing C API with zero tests is a release blocker. At
minimum: a single testthat file that uses `useDynLib(crs,
crs_nomad_solve)` + a small C/Rcpp shim, covering (a) a happy-path C
callback, (b) the R callback bridge, (c) every validation rejection
branch, (d) the option-mirror logic, (e) multi-start with both
caller-provided and crs-generated starts, (f) `read_nomad_opt_file =
1`, (g) bound coercion for integer/binary/categorical inputs, (h)
counter reporting (`callback_evaluations` vs failures).

Without this suite, every defect above had no chance of being caught
before publication.

---

## S3 — Medium-severity defects

### S3-1. Silent categorical-to-integer degradation

`map_input_type` maps `CRS_NOMAD_INPUT_CATEGORICAL` to NOMAD's `'I'`
(integer) with no diagnostic (snomadr_native.cpp:138–143).
`snomadr.cpp:225–226` emits `Rprintf("Warning: categorical input type
is mapped to integer in embedded NOMAD4.\n")` the first time. The
native API exposes a `CRS_NOMAD_INPUT_CATEGORICAL` enum value, so
callers will naturally assume it does something distinct.

Fix: either retire the enum value (and update the header to document
what's actually supported), or emit a one-shot diagnostic when
`quiet == 0`. A status field bit ("categorical_downgraded") in
`crs_nomad_result` would be the cleanest signal.

### S3-2. `random_seed == 0` has overloaded semantics

For start generation (`build_starting_points`):
`random_seed == 0` ⇒ use `std::time(nullptr)` (non-deterministic).

For NOMAD's SEED option (`apply_problem_parameters`):
`random_seed == 0` ⇒ no SEED set ⇒ NOMAD uses its built-in default
(typically 0, deterministic).

So when a caller wants non-deterministic behaviour by passing 0:
- if they let crs generate starts, generation is non-deterministic but
  NOMAD's internal search is deterministic;
- if they supply starts, NOMAD's internal search is deterministic
  (start_count > 1 path zeroes `random_seed` in `run_problem`, so
  even crs-generated starts result in default NOMAD SEED).

The header documents this for the "starts is NULL and start_count > 1"
case but does not warn about the asymmetry in other cases. End users
will get surprising determinism.

Fix: treat `random_seed == 0` as a single "engine chooses" sentinel
across both subsystems, or define an explicit
`CRS_NOMAD_SEED_NONDETERMINISTIC` / `_DEFAULT` sentinel. Document
clearly.

### S3-3. Code duplication with `snomadr.cpp`

`run_flag_message`, `suppress_native_nomad_display`,
`finite_or_default`, `size_to_int_saturated`, the input/output type
mapping switch statements, the `nomad.opt` file parser, the
`fill_random_start` body, the `build_starting_points` random+overwrite
pattern, and the `coerce_x0_value` body all duplicate code in
`snomadr.cpp`. They are *similar* but *not identical* (S2-4 is one
example of where they diverge).

Fix: factor into a shared `inst/include` (or `src/`) helper TU so both
paths are guaranteed to behave the same. The cost is one extra header
and one extra translation unit; the win is that S2-3, S2-4, S2-5,
S2-6, and S3-1 all collapse to single-source bugs.

### S3-4. `apply_native_options` heap-allocates per option

Line 526: `const std::string line = std::string(option->name) + " " +
option->value;`. Two allocations per option × number of options. Not a
hot path, but trivial to fix by using `addNomadStringParam(pb,
option->name, option->value)` first and falling back to
`addNomadParam(joined)` only when the typed setter is unsuitable —
which is what snomadr.cpp already does at lines 462–465 and 495–500.

### S3-5. `iterations` field never populated

`result->iterations = 0` in `initialize_result` (line 65); never
written elsewhere. snomadr.cpp does the same (line 949) but at least
the legacy contract doesn't promise iterations. The native header
advertises an `iterations` field; consumers will assume it carries
real information. Either populate from NOMAD's algorithm counters or
drop the field (since `struct_size` and reserved slots let you grow
the struct later).

### S3-6. `solve_final_problem` regenerates random starts even when
`start_count == 1` and `starts == NULL`

Line 620 short-circuits `start_count <= 0` to the single-start path,
but `start_count == 1` with `starts == NULL` runs through the random
generator (and then overwrites slot 0 with x0). One wasted RNG draw
per dimension. Trivial perf wart but reveals confusion in the
intended control flow.

Fix: short-circuit `start_count <= 1` to the single-start path.

### S3-7. `is_bad_number` allows ±Inf in x0

`is_bad_number(x)` returns `std::isnan(x)` only (line 196). Bounds
default to ±Inf so accepting Inf there is correct, but `x0[i] = +Inf`
will pass validation and reach `coerce_x0_value`, which only clamps
to upper if upper is finite. If upper is also Inf, the starting point
is Inf — NOMAD will likely reject it as run flag -3 ("initial point
failed to evaluate").

Fix: in validation, reject `!std::isfinite(x0[i])` (after permitting
bound-only Infs).

### S3-8. `validate_problem` does not check `bb_input_type` and
`bb_output_type` element ranges against the enum

It checks individual values via `map_input_type` /
`map_output_type` and rejects unknown codes, which is fine — but the
loop runs `n + m` switches. If the caller's array is null-deref'd
because `n` is huge but the array short (caller bug), validation
walks off into UB. The header treats those fields as borrowed
buffers; validation has to trust the caller's `n` and `m`. Acceptable
for a C ABI, but worth noting in the header's ownership prose.

---

## S4 — Low / polish

### S4-1. `solve_final_problem` naming is vestigial

Mentions "final" with no contrasting "initial". Holdover from v1/v2
days. Rename to `run_solve` or just inline.

### S4-2. Redundant `std::max(1, problem->start_count)` calls

Lines 628, 630, 640 all clamp; pick one local and reuse.

### S4-3. `result->blackbox_evaluations` is just a copy of
`result->callback_evaluations`

Lines 699–700 assign identical values. Either drop one field or have
`blackbox_evaluations` mean something distinct (e.g. NOMAD-reported
evaluations minus cache hits).

### S4-4. `initialize_result` writes NaN into solution/outputs *before*
validation

Lines 69–78. If `solution_len` is bogus (caller bug, e.g. a wildly
large value), this writes out of bounds before `validate_problem`
rejects the call. The header does require the caller to set
`solution_len` correctly, so this is contract-permissible — but
defensive code would gate the init on `solution_len >= 0 &&
solution_len <= reasonable_cap` or do the init after validation.

### S4-5. `solve_final_problem_with_starts` could be ~30 lines shorter
with RAII wrappers around `NomadProblem` and `NomadResult`

Mentioned in S1-5; also a structural improvement worth doing for
readability.

### S4-6. `crs_nomad_option`'s value is always a string

snomadr.cpp accepts both numeric and string option vectors and dispatches
to `addNomadValParam` / `addNomadDoubleParam` / `addNomadStringParam`.
The native API forces stringification at the call site, which means a
package author with an `int` MAX_BB_EVAL has to format it. Consider a
tagged-value struct (`int`/`double`/`string` discriminator) or a
parallel `crs_nomad_int_option` / `crs_nomad_double_option` set. Not
critical — string-only is a clean ABI — but worth deciding deliberately.

### S4-7. Header documents that result memory is caller-owned but the
ownership of `crs_nomad_r_callback` is implicit

Spell out: `crs_nomad_r_callback` is borrowed for the call's
duration; crs may dereference and follow its SEXP pointers; crs does
not retain copies.

### S4-8. `run_flag_message` doesn't surface NOMAD's own error text

When NOMAD fails parameter validation (`run_flag == -7`), users get
"NOMAD rejected the supplied parameters" with no detail. NOMAD's C++
API has a way to fetch the last error string. Plumbing it through
`crs_nomad_result.message` would materially help debugging.

### S4-9. `CRS_NOMAD_OK == 0` collision with snomadr's run-flag 0

Native's `result->status == CRS_NOMAD_OK` (0) means success. NOMAD's
own run flag 0 means "budget spent with feasible point" — also a
success but a *different* success. The code translates correctly but
the dual meaning of 0 across `status` and `nomad_run_flag` is a
documentation hazard.

### S4-10. Missing `extern "C"` cleanup in the `.cpp`

`snomadr_native.cpp:780` is `extern "C" int crs_nomad_solve(...)`. ✓
Header has `#ifdef __cplusplus extern "C"`. ✓
However, `#ifdef length`/`#ifdef error`/`#ifdef match` at the top of
the cpp guards against R.h's name-clashing macros — fine — but those
macros are typically `#undef`-guarded in `<Rinternals.h>`. Verify
they're actually defined when this TU is compiled in the current
toolchain; if not, the `#undef`s are dead code that should be removed.

---

## Performance and leanness

The hot path is `native_eval_single → native_eval_c` /
`native_eval_r`, called once per blackbox evaluation. The C path is
already lean: a single function-pointer call, an O(m) finite-check
loop. Two cheap wins remain:

(P-1) `native_eval_single` does `++context->callback_evaluations`
unconditionally; with S2-3's fix, restructure so the increment is at
the end alongside `*count_eval = true` (one branch becomes one fused
write).

(P-2) The finite-check loop at lines 364–369 runs after every callback.
For an EB constraint that legitimately returns +Inf to signal
infeasibility, treating it as a callback failure is too strict —
snomadr does the same finite check via NOMAD's own machinery (NOMAD
treats Inf as infeasible without calling it a failure). Current code
will report `callback_failures++` for any +Inf output, which
combined with S2-3 will overcount failures and then convert the run
to `CRS_NOMAD_CALLBACK_FAILURE` even though the solver may have
succeeded. **This is closer to a correctness defect — promote to
S2-8.** Fix: let NOMAD handle Inf via its EB/PB output classifications;
only reject NaN.

The R path is dominated by `Rf_allocVector` + `R_tryEval` per call.
Negligible improvements available; the cost is intrinsic to the
R-bridge contract.

Outside the hot loop, allocations per solve are minor:
`input_type_string`, `output_type_string`, two `std::vector<double>`
for best_x/best_out, and the starts vector. All scale linearly in
problem dimension. No micro-optimisations worth the readability cost.

---

## Engineering structure

Header design is good: versioned, growable, opaque, callback modes
explicit. The single registered callable is correct policy.

Implementation:

- Single 815-line cpp with one anonymous namespace and one extern "C"
  entry. Reasonable for the surface area but would benefit from
  splitting into `snomadr_native.cpp` (entry, validation, orchestration),
  `snomadr_native_eval.cpp` (C/R callback bridges, finite checks), and
  a shared `nomad_helpers.cpp` with snomadr (S3-3).
- Error reporting is uniform through `fail()` — good.
- Magic numbers for NOMAD run flags would benefit from named constants
  matching NOMAD's internal enum, declared once at the top of the TU.
- Naming is consistent (snake_case, `crs_nomad_*`).
- No dead code in the TU, but `solve_final_problem` reads as if there
  were a non-final variant — vestigial (S4-1).

## Cross-check against `snomadr()` behavioural parity

Differences observed (severity in brackets):

1. `callback_evaluations` overcount on failures [S2-3].
2. No `R_CheckUserInterrupt` [S2-1].
3. No categorical-to-integer warning [S3-1].
4. `nomad.opt` parsing doesn't trim leading whitespace [S2-4].
5. `MAX_BB_EVAL` mirror requires strict decimal integer string [S2-5].
6. Option-name matching is case-sensitive [S2-6].
7. `random_seed` semantics differ between start generation and NOMAD
   SEED [S3-2].
8. Native validates "at least one OBJ output present"; snomadr defers
   to NOMAD. Stricter is fine; document it.
9. `INF` outputs become callback failures [S2-8 (promoted)].

Otherwise: starting-point generation algorithm (random-fill then
overwrite slot 0 with coerced x0) matches; bound coercion matches;
result loading matches; cache/total evaluation fallbacks match.

---

## Recommended next actions, ranked

1. Build the unit test suite for `crs_nomad_solve` (S2-7). Without
   it, every fix below is undefended.
2. Fix the exception/longjmp safety hole: `R_PreserveObject` for the R
   callback SEXPs, RAII for NOMAD handles, busy-flag release in a
   cleanup-guaranteed path (S1-1, S1-2, S1-5).
3. Remove the +Inf-as-failure trap (S2-8).
4. Add `R_CheckUserInterrupt()` (S2-1).
5. Harden option matching and the nomad.opt parser (S2-4, S2-5, S2-6).
6. Fix callback counter accounting (S2-3).
7. Factor shared helpers between `snomadr.cpp` and
   `snomadr_native.cpp` (S3-3) — this is where the long-term defect
   reduction lives.
8. Coordinate process-wide state between the two NOMAD entry points
   (S1-4).
9. Document or remove `crs_nomad_result.iterations`,
   `blackbox_evaluations` aliasing, categorical handling, and seed
   semantics (S3-1, S3-2, S3-5, S4-3, S4-9).
10. Add an example consumer (a tiny demo package or a vignette code
    snippet) so external maintainers can verify the contract.

---

## Addendum — findings derived from the NOMAD 4 user guide

Source: NOMAD 4 user guide,
<https://nomad-4-user-guide.readthedocs.io/en/latest/>. Per-finding
links to the relevant section are inline below; the full set of pages
consulted is gathered in the *Sources* block at the end of this
addendum.

The user guide surfaces six material defects in the native API that the
source review alone could not see, plus several confirmations of
earlier findings.

### S1-6 (new — critical). `NB_THREADS_PARALLEL_EVAL` > 1 in R-callback mode is catastrophic UB

The [*Parallel evaluations*](https://nomad-4-user-guide.readthedocs.io/en/latest/AdvancedFunctionalities.html#parallel-evaluations)
page states explicitly:

> "When OpenMP is available, the user MUST provide the number of
> threads for the parallel evaluations `NB_THREADS_PARALLEL_EVAL` to
> efficiently access the computer cores. If this parameter is not set,
> OpenMP uses a single thread. The evaluations of trial points stored
> in a queue are dispatched to these threads."

The native API passes raw user-supplied options through `addNomadParam`
without filtering. A caller in `CRS_NOMAD_CALLBACK_R` mode who sets
`NB_THREADS_PARALLEL_EVAL 8` (perfectly reasonable for a heavy
objective) will have NOMAD dispatch the evaluator onto N OpenMP worker
threads, each of which then calls `native_eval_r`, each of which then
calls the R API. Concurrent R API calls are unambiguous undefined
behaviour — instant crash, corrupted SEXPs, or both.

This is worse than the thread-safety concern flagged in S1-3 because
it requires no caller error to trigger: the user is following NOMAD's
own documentation when they set the parameter.

Fix: in `apply_native_options` / `apply_problem_parameters`, when
`callback_mode == CRS_NOMAD_CALLBACK_R`:

1. Force `NB_THREADS_PARALLEL_EVAL 1` after applying caller options, or
2. Reject the call with `CRS_NOMAD_INVALID_INPUT` and a clear message
   if the caller asked for >1 threads.

Option (1) is friendlier; option (2) is louder. Either is required.

The C-callback path is fine in principle but the caller's blackbox
must itself be thread-safe — document this.

### S2-9 (new). The API offers no way to request non-deterministic NOMAD search

The user guide's [`SEED`](https://nomad-4-user-guide.readthedocs.io/en/latest/AdvancedFunctionalities.html#seed)
page:

> "The user can change the sequence of directions by setting `SEED` to
> a positive integer or `-1`. If `-1` or `DIFF` is entered the seed is
> different for each run (PID is used)."

The native header types `random_seed` as `unsigned int`. A caller
*cannot* express `-1`. The implementation does not pass `SEED` when
`random_seed == 0` (snomadr_native.cpp:588), so the field has only one
useful state: "set NOMAD to a specific positive value." There is no
way through the structured surface to request PID-seeded behaviour.
Callers can still pass `"SEED -1"` through the options array, but the
header advertises `random_seed` as the way to control seeding —
silently.

This compounds S3-2 ("`random_seed == 0` has overloaded semantics"):
the field is now revealed to have *zero* useful semantics for the
non-deterministic case.

Fix: define a sentinel — either reinterpret `random_seed == 0` as
"NOMAD chooses (default deterministic 0)" and add a `quiet_seed` /
`use_pid_seed` / similar flag, or change `random_seed` to a signed
type with `-1` meaning "PID seed" matching NOMAD's surface. The latter
is a breaking ABI change, so for v1 either add a new field through
one of the reserved slots, or document that callers wanting
non-determinism must pass `"SEED -1"` via options.

### S2-10 (promoted from S3-1). NOMAD 4 does not support categorical variables at all

The user guide's [`BB_INPUT_TYPE`](https://nomad-4-user-guide.readthedocs.io/en/latest/HowToUseNomad.html#bb-input-type)
page:

> "Input types `t` are values in `R, B, I`."
>
> Note: "Categorical variables are not yet supported in NOMAD 4 but
> are available in NOMAD 3. Some work by a PhD student is currently
> being done to reintroduce this feature in NOMAD 4."

The native header still exposes `CRS_NOMAD_INPUT_CATEGORICAL` as a
first-class enum value. The implementation silently downcasts it to
`'I'` (integer). This is an active misrepresentation of capability,
made worse by being a brand-new public API — a downstream package
author will reasonably assume an exposed enum value does something
distinct.

snomadr keeps the same downcast but at least prints a warning. The
native API does not.

Fix: remove `CRS_NOMAD_INPUT_CATEGORICAL` from the enum (it's API v1,
the only break is symbolic; the wire value is currently 2 — repurpose
or reserve), or document it in the header as
"deprecated/reserved-pending-NOMAD-support" and emit a structured
diagnostic via the result's `message` / a new reserved status bit. Do
not silently downcast.

### S3-9 (new). NOMAD's cache is a process-wide singleton

The user guide's [`CACHE_FILE`](https://nomad-4-user-guide.readthedocs.io/en/latest/HowToUseNomad.html#cache-file)
and library-mode [*Access to solution and optimization data*](https://nomad-4-user-guide.readthedocs.io/en/latest/LibraryMode.html#access-to-solution-and-optimization-data)
sections describe `NOMAD::CacheBase` as a singleton accessed via
`getInstance()`. The C interface allocates
fresh `NomadProblem` and `NomadResult` instances per
`solveNomadProblem` call, but the underlying cache singleton persists
across calls in the same process.

Consequence: two sequential calls to `crs_nomad_solve` from the same R
session will share evaluated points unless one calls back to NOMAD's
cache management. If the second call has a different blackbox or
bounds, NOMAD may return cached evaluations from the wrong problem
(it keys on point coordinates, not problem identity).

This affects `snomadr()` identically and was probably never the
intended user-facing semantics of a stateless C ABI.

Fix: either set `EVAL_USE_CACHE false` by default (caller can opt back
in), or clear the singleton between calls — NomadStdCInterface
doesn't expose a clear hook directly, but `MainStep::resetCache()`
or calling `CacheBase::getInstance()->clear()` from the C++ side
during teardown should work. Document either way.

### S3-10 (new). Missing output types relative to NOMAD's surface

Per [`BB_OUTPUT_TYPE`](https://nomad-4-user-guide.readthedocs.io/en/latest/HowToUseNomad.html#bb-output-type),
NOMAD 4 supports `OBJ`, `PB`/`CSTR`, `EB`, `F`, `CNT_EVAL`, and
`NOTHING`/`EXTRA_O`/`-` as `BB_OUTPUT_TYPE` values. The native API
exposes only `OBJ`, `PB`, `EB`. Missing:

- `F` (filter constraint) — different mathematical treatment from PB.
- `CNT_EVAL` — lets NOMAD decide whether a row counts as an evaluation.
- `NOTHING` / `EXTRA_O` / `-` — diagnostic/passthrough outputs that
  appear in cache and history but don't affect the algorithm.

snomadr.cpp also only handles `OBJ`/`PB`/`EB`, but the R-side
`bbout` argument is loose enough that a sophisticated caller could
work around. The native API is a fresh ABI with no such escape hatch.

Fix: extend `crs_nomad_output_type` to cover at least `F` and
`NOTHING`/`EXTRA_O`; add the wire-level mapping in `map_output_type`.
Use reserved enum slots if v1 ABI freeze is in effect.

### S3-11 (new). Block-evaluation callback is not surfaced

NomadStdCInterface defines `Callback_BB_block` for evaluating
multiple trial points per call (advertised under `BB_MAX_BLOCK_SIZE`
in the [*Blackbox evaluation of a block of trial points*](https://nomad-4-user-guide.readthedocs.io/en/latest/AdvancedFunctionalities.html#blackbox-evaluation-of-a-block-of-trial-points)
section). The native API only wires `bb_single`
(snomadr_native.cpp:673). Users who set `BB_MAX_BLOCK_SIZE > 1` will
have NOMAD coerce back to single, silently losing the per-call
amortisation. For parallel-amenable blackboxes (the typical reason
to use block eval) this is a substantial performance loss.

Fix: add a `crs_nomad_eval_block_fn` typedef and a
`callback_mode == CRS_NOMAD_CALLBACK_C_BLOCK` value; wire through
`createNomadProblem(nullptr, native_eval_block, n, m)` when set.
Trivial extension once the design slot is reserved.

### S4-11 (new). `nomad.opt` parser doesn't strip inline comments

The [*Basic parameters description*](https://nomad-4-user-guide.readthedocs.io/en/latest/HowToUseNomad.html#basic-parameters-description)
section:

> "The remaining content of a line is ignored after the character
> `#`."

`apply_nomad_opt_file` only recognises `#` at the first non-whitespace
position (snomadr_native.cpp:550) and passes the full line otherwise.
NOMAD's own `addNomadParam` parser likely handles inline `#`
itself — verify in NomadStdCInterface source — but to match snomadr
and the documented file format, strip inline comments client-side
before forwarding.

### S4-12 (new). Bounds-as-"undefined" convention not documented in header

The [*Basic parameters description*](https://nomad-4-user-guide.readthedocs.io/en/latest/HowToUseNomad.html#basic-parameters-description)
section (see also [`LOWER_BOUND` and `UPPER_BOUND`](https://nomad-4-user-guide.readthedocs.io/en/latest/HowToUseNomad.html#lower-bound-and-upper-bound)):

> "The strings `-`, `inf`, `-inf` or `+inf` are accepted to enter
> undefined real values (NOMAD considers ±∞ as an undefined value)."

The native API forces `lower` and `upper` to be non-null pointers and
documents nothing about how to signal "no bound." A caller using
`-DBL_MAX` instead of `-INFINITY` will get NOMAD treating that as a
real bound, with bizarre numerical consequences.

Fix: header prose — "set element to `-INFINITY` / `INFINITY` to signal
unbounded; ±∞ is treated by NOMAD as 'undefined'." One sentence,
fixes a class of caller errors.

### S4-13 (new). No convenience surface for `MAX_TIME` or `LH_SEARCH`

NOMAD's [basic algorithmic parameter table](https://nomad-4-user-guide.readthedocs.io/en/latest/HowToUseNomad.html#algorithmic-parameters)
lists `MAX_TIME` and `LH_SEARCH` (Latin Hypercube initial search)
alongside `MAX_BB_EVAL` as headline stopping/initialisation controls. The native API hoists
`max_eval` to a struct field but leaves the other two to the option
array. Asymmetric — and `MAX_TIME` is the standard escape hatch when
a user cannot estimate `MAX_BB_EVAL`. Consider adding `int
max_time_seconds` and `int lh_initial_eval` to the next struct
revision (reserved slots are available).

### Confirmations of earlier findings

The guide also tightens or confirms:

- **S2-6** (case-sensitive option matching): NOMAD parameter names are
  documented as case-insensitive. Native API's `strcmp` for
  `MAX_EVAL`/`MAX_BB_EVAL` is therefore inconsistent with the
  backend.

- **S2-1** (no Ctrl-C interrupt) is more severe than first noted: when
  combined with `MAX_TIME` not being a struct convenience field
  (S4-13) and the parallel-eval hazard (S1-6), a user running with a
  long-running objective has no clean way to abort.

- **S2-8** (Inf-as-failure) is reinforced by the guide's description of
  `PB`/`EB` constraint handling — NOMAD genuinely expects ±∞ from
  the blackbox as a signal of infeasibility, and converting that to
  a callback failure is wrong by design.

- **S1-3** / **S1-1** (R API outside main thread / longjmp safety) are
  reinforced by NOMAD's explicit warning in the
  [library-mode *Setting parameters*](https://nomad-4-user-guide.readthedocs.io/en/latest/LibraryMode.html#setting-parameters)
  section that "NOMAD routines may throw `C++` exceptions, so it is
  recommended that you put your code into a `try` block." The native
  API has the outer try/catch but not the inner-frame leak
  protection (S1-5).

### Updated action ranking

Inserting the new findings into the ranked list:

1. Build the unit test suite (S2-7) — still first.
2. Force / reject `NB_THREADS_PARALLEL_EVAL > 1` in R-callback mode
   (S1-6 — new). This is the loudest crash vector and trivial to
   patch.
3. Fix the exception/longjmp safety hole (S1-1, S1-2, S1-5).
4. Remove `CRS_NOMAD_INPUT_CATEGORICAL` or mark it reserved (S2-10).
5. Remove the +Inf-as-failure trap (S2-8).
6. Add `R_CheckUserInterrupt()` (S2-1).
7. Harden option matching and the nomad.opt parser (S2-4, S2-5,
   S2-6, S4-11).
8. Fix callback counter accounting (S2-3).
9. Add a non-deterministic seed mechanism (S2-9).
10. Reset / document the NOMAD cache singleton (S3-9).
11. Factor shared helpers with snomadr.cpp (S3-3).
12. Surface missing output types and block eval (S3-10, S3-11).
13. Coordinate process-wide state between the two NOMAD entry points
    (S1-4).
14. Document or expand the convenience surface — bounds-as-undefined
    (S4-12), `MAX_TIME` / `LH_SEARCH` (S4-13), seed semantics
    (S3-2), iterations / aliasing (S3-5, S4-3), categorical
    handling and `OK == 0` overload (S4-9).

---

### Sources

NOMAD 4 user guide pages consulted for this addendum:

- [NOMAD 4 User Guide — home](https://nomad-4-user-guide.readthedocs.io/en/latest/)
- [NOMAD usage — basic parameters](https://nomad-4-user-guide.readthedocs.io/en/latest/HowToUseNomad.html)
  - [`BB_INPUT_TYPE`](https://nomad-4-user-guide.readthedocs.io/en/latest/HowToUseNomad.html#bb-input-type)
  - [`BB_OUTPUT_TYPE`](https://nomad-4-user-guide.readthedocs.io/en/latest/HowToUseNomad.html#bb-output-type)
  - [`LOWER_BOUND` and `UPPER_BOUND`](https://nomad-4-user-guide.readthedocs.io/en/latest/HowToUseNomad.html#lower-bound-and-upper-bound)
  - [Algorithmic parameters (`MAX_BB_EVAL` / `MAX_TIME` / `LH_SEARCH`)](https://nomad-4-user-guide.readthedocs.io/en/latest/HowToUseNomad.html#algorithmic-parameters)
  - [`CACHE_FILE`](https://nomad-4-user-guide.readthedocs.io/en/latest/HowToUseNomad.html#cache-file)
- [Optimization in library mode](https://nomad-4-user-guide.readthedocs.io/en/latest/LibraryMode.html)
  - [Setting parameters (C++ exception warning)](https://nomad-4-user-guide.readthedocs.io/en/latest/LibraryMode.html#setting-parameters)
  - [Access to solution and optimization data](https://nomad-4-user-guide.readthedocs.io/en/latest/LibraryMode.html#access-to-solution-and-optimization-data)
  - [C interface](https://nomad-4-user-guide.readthedocs.io/en/latest/LibraryMode.html#c-interface)
- [Advanced functionalities](https://nomad-4-user-guide.readthedocs.io/en/latest/AdvancedFunctionalities.html)
  - [`SEED`](https://nomad-4-user-guide.readthedocs.io/en/latest/AdvancedFunctionalities.html#seed)
  - [`EVAL_OPPORTUNISTIC`](https://nomad-4-user-guide.readthedocs.io/en/latest/AdvancedFunctionalities.html#eval-opportunistic)
  - [Blackbox evaluation of a block of trial points](https://nomad-4-user-guide.readthedocs.io/en/latest/AdvancedFunctionalities.html#blackbox-evaluation-of-a-block-of-trial-points)
  - [Parallel evaluations](https://nomad-4-user-guide.readthedocs.io/en/latest/AdvancedFunctionalities.html#parallel-evaluations)
- [Complete list of parameters (appendix)](https://nomad-4-user-guide.readthedocs.io/en/latest/Appendix.html#complete-list-of-parameters)

---

End of review.
