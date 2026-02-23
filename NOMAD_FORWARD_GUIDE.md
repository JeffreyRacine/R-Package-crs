# NOMAD Forward Guide for `crs`

Last updated: 2026-02-23

This document is the operational guide for NOMAD work in `crs`.

Companion technical map:

- `/Users/jracine/Development/crs/NOMAD_INTEGRATION_MAP.md`

## Current baseline (must preserve)

- Embedded NOMAD version wired into `crs`: `4.5.0`
- Embedded source trees:
  - `/Users/jracine/Development/crs/src/nomad4_src`
- Legacy retained source trees (no longer wired into build):
  - `/Users/jracine/Development/crs/src/nomad_src`
  - `/Users/jracine/Development/crs/src/sgtelib_src`
- R entry wrapper:
  - `/Users/jracine/Development/crs/R/snomadr.R`
- Native bridge:
  - `/Users/jracine/Development/crs/src/snomadr.cpp`
- Registration and DLL exposure:
  - `/Users/jracine/Development/crs/src/crs_init.c`
  - `/Users/jracine/Development/crs/NAMESPACE`

## Primary usage paths in package code

- `crs(..., cv="nomad")` dispatches to:
  - `/Users/jracine/Development/crs/R/frscvNOMAD.R`
  - `/Users/jracine/Development/crs/R/krscvNOMAD.R`
- Additional NOMAD-dependent paths:
  - `/Users/jracine/Development/crs/R/clsd.R`
  - `/Users/jracine/Development/crs/R/np.regression.glp.R`

## Working rule for future edits

1. Treat NOMAD integration as a three-layer system:
   - R API layer (`snomadr.R` and callers)
   - C++ bridge layer (`snomadr.cpp/h`)
   - Vendored solver layer (`nomad_src`, `sgtelib_src`, `Makevars`)
2. Do not change one layer without checking interface compatibility against adjacent layers.
3. Keep numerical behavior stable in fit/predict/CV pathways unless intentionally changed and documented.

## Recommended modernization sequence

1. Lock a baseline
   - Build and run current package checks.
   - Record representative NOMAD-driven outputs from:
     - `frscvNOMAD`
     - `krscvNOMAD`
     - `glpcvNOMAD`
2. Isolate bridge contract
   - Document exact expected structure passed from R to `.Call(...)`:
     - `eval.f`, `n`, `bbin`, `bbout`, `x0`, bounds, `options`, `nmulti`, seed.
   - Confirm all bridge return fields consumed by R (`status`, `message`, `bbe`, `iterations`, `objective`, `solution`).
3. Prepare staged upgrade branch for solver source
   - Introduce newer NOMAD sources in a separate staging location first.
   - Adapt `Makevars` object list and include flags only after bridge compile targets are clear.
4. Port bridge API deliberately
   - Reconcile API differences for `Parameters`, solver run lifecycle, callbacks, and status codes.
   - Keep R-facing return contract unchanged initially.
5. Validate parity before optimization
   - Run fixed-seed and varying-seed comparisons for main NOMAD paths.
   - Confirm no unintended regressions in objective values and selected hyperparameters.

## Minimum acceptance checks for any NOMAD change

1. Build succeeds from clean objects.
2. `R CMD check` has no new warnings/notes attributable to NOMAD edits.
3. NOMAD paths run end-to-end:
   - `frscvNOMAD`
   - `krscvNOMAD`
   - at least one `glpcvNOMAD` call
4. Return-object compatibility remains intact for existing R callers.
5. Numerical drift is either negligible or explicitly justified.

## Risks to watch

- Option handling mismatch between R and C++ (`get.option.types.R` and parameter parsing)
- Status-code or message mapping drift in `snomadr.h` / `snomadr.cpp`
- Multi-start behavior changes (`smultinomadRSolve`) affecting CV stability
- Build graph breakage due to incomplete object list updates in `src/Makevars`

## Immediate next action (recommended)

- Add a small reproducibility harness dedicated to NOMAD paths (single script) and store baseline outputs before any version upgrade work begins.

## Roadmap update (2026-02-23)

The cross-package dependency concern around `crs`/`np` for NOMAD is no longer treated as a blocker under the current roadmap.

Current roadmap intent:

1. `crs::npglpreg` is expected to be deprecated once the migrated/extended `np::npreg` pathway is the canonical replacement.
2. The replacement functionality is delivered through `npreg` using:
   - `regtype = "lp"`
   - `basis = ...`
   - `degree = ...`
3. This migrated implementation has recently been validated to produce numerically identical results relative to the prior `npglpreg` behavior.

Decision point now framed for future work:

- Determine whether to augment `npreg` to use the `crs` NOMAD solver path to deliver the `npglpreg`-style optimization workflow currently housed in `crs`.

Implication for NOMAD planning:

- Near-term NOMAD work should preserve `crs` solver behavior and keep the bridge stable so it can be reused if the `npreg` augmentation path is selected.

## Suggested practical sequence: upgrade NOMAD in `crs` with strict parity/performance gates

This sequence is designed to enforce the same discipline used in your modernization workflow (explicit pre/post, fixed+varying seeds, mean/median deltas, objective parity, `/tmp` artifacts).

### Phase 0: Baseline freeze (no NOMAD changes yet)

1. Create a dedicated branch and capture baseline commit hash.
2. Build/check current package from tarball:
   - `cd /Users/jracine/Development`
   - `R CMD build crs`
   - `R CMD check --as-cran crs_<version>.tar.gz`
3. Save logs:
   - `/tmp/crs_nomad_baseline_build.log`
   - `/tmp/crs_nomad_baseline_check.log`

Gate:

- Baseline check status is known and archived before any NOMAD edits.

### Phase 1: Build reproducible NOMAD harnesses

Create two harness scripts in `crs/benchmarks/nomad/` (or equivalent):

1. `run_nomad_smoke.R` (small `n`, small repeats)
   - Covers:
     - `frscvNOMAD`
     - `krscvNOMAD`
     - `npglpreg`/`glpcvNOMAD` path
   - Emits:
     - objective values
     - selected degrees/segments/lambda/bandwidths
     - elapsed times
2. `run_nomad_bench.R` (larger runs for stable timing)
   - Same outputs with replicated runs for mean/median summaries.

Required run policies:

- Fixed-seed run set.
- Varying-seed run set.
- Explicitly report mode as `serial`.

Gate:

- Smoke harness passes on baseline and produces structured outputs.

### Phase 2: Pre-change measurements (reference dataset)

Run both harnesses before touching NOMAD source:

- `/tmp/crs_nomad_pre_smoke_raw.csv`
- `/tmp/crs_nomad_pre_smoke_summary.csv`
- `/tmp/crs_nomad_pre_bench_raw.csv`
- `/tmp/crs_nomad_pre_bench_summary.csv`
- `/tmp/crs_nomad_pre_parity.rds`

Gate:

- Pre-change metrics are complete (fixed+varying seed).

### Phase 3: NOMAD upgrade in isolation

1. Stage updated NOMAD source in a temporary directory under `crs/src` (do not overwrite blind).
2. Port `snomadr.cpp`/`snomadr.h` bridge API against new NOMAD API.
3. Update `src/Makevars` object/include graph.
4. Keep R-level `snomadr` return contract unchanged:
   - `status`, `message`, `bbe`, `iterations`, `objective`, `solution`.

Gate:

- Package compiles and loads; NOMAD smoke calls execute end-to-end.

### Phase 4: Post-change measurements and strict parity checks

Run identical harnesses with identical arguments:

- `/tmp/crs_nomad_post_smoke_raw.csv`
- `/tmp/crs_nomad_post_smoke_summary.csv`
- `/tmp/crs_nomad_post_bench_raw.csv`
- `/tmp/crs_nomad_post_bench_summary.csv`
- `/tmp/crs_nomad_post_parity.rds`

Compute comparisons:

- `/tmp/crs_nomad_prepost_perf_summary.csv`
  - mean percent change
  - median percent change
- `/tmp/crs_nomad_prepost_parity_summary.csv`
  - objective parity (`all.equal`/max abs diff)
  - selected parameter parity (degree/segments/lambda/bandwidth)
  - solution-vector deltas

Gate (strict):

1. No unexplained objective regressions.
2. Parameter/solution drift either:
   - negligible under tolerance, or
   - explicitly justified (e.g., solver ordering/API behavior change).
3. Performance impact quantified with mean+median deltas (fixed+varying seeds).

### Phase 5: Package QA and regression gates

1. Run package tests:
   - `R -q -e "devtools::test('/Users/jracine/Development/crs')"` (or equivalent)
2. Tarball-first check again:
   - `cd /Users/jracine/Development`
   - `R CMD build crs`
   - `R CMD check --as-cran crs_<version>.tar.gz`
3. Save logs:
   - `/tmp/crs_nomad_post_build.log`
   - `/tmp/crs_nomad_post_check.log`
   - `/tmp/crs_nomad_testthat.log`

Gate:

- No new unexplained warnings/notes/errors attributable to NOMAD changes.

### Phase 6: Decision checkpoint for roadmap alignment

If gates pass, choose one of:

1. Keep upgraded NOMAD fully internal to `crs` (short-term stability path).
2. Start extraction of solver bridge into standalone package with the same validated contract.
3. Prepare targeted `npreg` augmentation to consume equivalent solver behavior for `regtype="lp"` workflows.

Deliverable for checkpoint:

- Single markdown report (stored in repo) containing:
  - script names + exact invocations
  - fixed/varying seed policy
  - mean/median % deltas
  - parity results and tolerances
  - `/tmp` artifact inventory

## Checkpoint Addendum (2026-02-23, NOMAD4 bridge hardening + gates)

### Scope completed

1. `snomadr.cpp` bridge moved to and hardened for NOMAD4 C interface behavior.
2. `Makevars` wiring updated to NOMAD4 source graph.
3. Strict pre/post benchmark harnesses executed with isolated pre/post libraries.

### Notable bridge compatibility updates

1. Array-option handling was made dimension-aware and NOMAD4-safe.
2. Legacy poll options are translated:
   - `MIN_POLL_SIZE` -> `MIN_FRAME_SIZE`
   - `INITIAL_POLL_SIZE` -> `INITIAL_FRAME_SIZE`
3. Relative (`r...`) scalar values are parsed for array options.
4. Integer/granular variables enforce mesh/frame lower bounds needed by NOMAD4.
5. `MAX_EVAL` is auto-set from `MAX_BB_EVAL` when missing (legacy budget intent).
6. `EPSILON` values are sanitized against NOMAD4 precision lower bound.

### Gate artifacts (latest run set)

Pre baseline:

- `/tmp/crs_nomad_pre_20260223_104449_raw.csv`
- `/tmp/crs_nomad_pre_20260223_104449_summary.csv`

Post (latest candidate):

- `/tmp/crs_nomad_post_20260223_104449_v2_raw.csv`
- `/tmp/crs_nomad_post_20260223_104449_v2_summary.csv`

Comparisons:

- `/tmp/crs_nomad_perf_20260223_104449_v2.csv`
- `/tmp/crs_nomad_parity_20260223_104449_v2.csv`

Install/bench logs:

- `/tmp/crs_nomad_gate_20260223_104449.log`
- `/tmp/crs_nomad_gate_20260223_104449_v2.log`

### Gate summary (latest)

1. Performance improved materially versus earlier NOMAD4 attempt after `MAX_EVAL` compatibility capping.
2. Strict performance parity is still not met in several paths (`frscvNOMAD`, `krscvNOMAD`, `snomadr_basic_lib` fixed-seed set).
3. Objective parity is strong for `frscvNOMAD`; mixed/partial for other paths.
4. Parameter parity still differs in `npglpreg` relative to NOMAD3 baseline.

### Current conclusion

The NOMAD4 upgrade is wired and operational with hardened bridge compatibility, but final strict parity/performance acceptance gates are not yet fully satisfied and require additional tuning decisions.
