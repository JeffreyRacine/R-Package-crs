# NOMAD Forward Guide for `crs`

Last updated: 2026-02-23

This document is the operational playbook for future NOMAD work in `crs`.

## Canonical doc split (no overlap)

1. Architecture and wiring (static):
   - `/Users/jracine/Development/crs/NOMAD_INTEGRATION_MAP.md`
2. Current gate status, artifacts, and open issues (dynamic):
   - `/Users/jracine/Development/crs/NOMAD_STATUS_AND_GATES.md`
3. Harness mechanics and commands:
   - `/Users/jracine/Development/crs/benchmarks/nomad/README.md`
4. Complete NOMAD4 option reference used by `crs`:
   - `/Users/jracine/Development/crs/inst/nomad/NOMAD_4_5_0_OPTIONS_REFERENCE.md`
5. Legacy-to-current option mapping:
   - `/Users/jracine/Development/crs/NOMAD_391_TO_450_OPTION_MAP.md`

## Current strategic baseline

1. `crs` is wired to embedded NOMAD `4.5.0`.
2. `crs` keeps a NOMAD4-only source layout (legacy NOMAD3 tree removed).
3. `snomadr()` applies a NOMAD3.9.1-like compatibility profile via NOMAD4 option names when those options are not user-specified.
4. Near-term goal remains: stabilize NOMAD4 under strict parity/performance gates for existing `crs` behavior.
5. `frscvNOMAD` and `krscvNOMAD` now apply path defaults (`QUAD_MODEL_SEARCH=no`, `EVAL_QUEUE_SORT=DIR_LAST_SUCCESS`) unless user overrides.

## Roadmap alignment (`crs` vs `np`/`npRmpi`)

1. `crs::npglpreg` is expected to be deprecated after full migration to `np::npreg` (`regtype="lp"`, `basis=`, `degree=`).
2. NOMAD work in `crs` should preserve solver behavior so that `npreg` can reuse it if needed.
3. Cross-package dependency concern is currently non-blocking; correctness and reproducibility gates are primary.

## Non-negotiable acceptance gates

1. Build/install from clean objects must succeed.
2. No new NOMAD-attributable `R CMD check --as-cran` regressions.
3. Pre/post comparisons must be run with:
   - fixed seeds
   - varying seeds
   - mean/median performance deltas
   - objective/parameter parity summaries
4. Any residual drift must be explicitly documented with likely source and practical impact.

## Practical next-step priorities

1. Close remaining NOMAD4 parity/performance gaps recorded in status doc.
2. Avoid global option toggles that improve one path but regress another; use path-aware tuning if needed.
3. Keep return contract stable (`status`, `message`, `bbe`, `iterations`, `objective`, `solution`) while tuning internals.

## NOMAD4 solver-mode guidance

1. NOMAD4 exposes multiple standalone algorithms (`CS_OPTIMIZATION`, `NM_OPTIMIZATION`, `QP_OPTIMIZATION`, `QUAD_MODEL_OPTIMIZATION`, etc.), but many are either specialized or not suitable for these mixed-variable CV objectives.
2. Practical recommendation for `crs`/`npglpreg`:
   - keep MADS as the primary optimization engine,
   - tune search-method options path-by-path,
   - avoid switching to standalone solver modes without explicit parity gates.
3. OpenMP-gated modes (`PSD_MADS_OPTIMIZATION`, `COOP_MADS_OPTIMIZATION`) are unavailable in the current `crs` build configuration and should not be considered default candidates.
4. Global-optimum practicality:
   - NOMAD is a local/direct-search framework; there is no strict global-optimum guarantee for these nonconvex mixed-variable CV objectives.
   - fastest reliable practice is path-tuned MADS plus controlled restarts (`nmulti`) and seed-policy checks rather than switching to standalone modes that degrade objective quality.

## Current default recommendation matrix

Based on the expanded solver-space sweeps documented in status:

1. `frscvNOMAD`
   - keep current MADS path defaults (`QUAD_MODEL_SEARCH=no`, `EVAL_QUEUE_SORT=DIR_LAST_SUCCESS`)
   - optional `EVAL_QUEUE_SORT=LEXICOGRAPHICAL` can be tested as an experimental profile (small mean speed gain, not consistently better in every cell)
2. `krscvNOMAD`
   - keep current MADS path defaults
   - do not switch to standalone solvers for default use (large speed gains came with objective-risk counts)
3. `npglpreg`
   - keep current MADS compatibility defaults
   - `QUAD_MODEL_SEARCH_BOX_FACTOR=1.0` and `QUAD_MODEL_BOX_FACTOR=1.0` (`np_quadbox1`) is the safest optional variant tested; speed impact is small/mixed

## What was tuned beyond solver toggles

The investigation covered, at minimum:

1. standalone solver families (`CS`, `NM`, `QP`, `RANDOM`)
2. MADS queue/search variants (`QUAD_MODEL_SEARCH`, `EVAL_QUEUE_SORT`, quad-box factors)
3. restart and budget controls (`nmulti`, `max.bb.eval`)
4. both fixed-seed and varying-seed comparisons on simple and harder mixed-variable scenarios

## Exception-only modes (important)

1. In this build, some modes are unavailable or require setup not present in `crs` test paths:
   - `VNS_MADS_OPTIMIZATION` (not implemented here)
   - DiscoMADS without revealing outputs
2. When such exception profiles are tested, isolate them from main comparisons (separate process/session) to avoid contaminated multi-profile runs.

## NOMAD release-update playbook

Use this sequence whenever a new NOMAD4.x release is integrated.

1. Vendor update and build wiring
   - update `src/nomad4_src` to the target upstream release
   - verify `src/Makevars` still compiles only from `src/nomad4_src`
   - keep clean-break policy: do not reintroduce `src/nomad_src`
2. Option contract policy
   - keep direct NOMAD4 option names only (`MIN_FRAME_SIZE`, `INITIAL_FRAME_SIZE`, etc.)
   - do not add legacy alias mapping (`MIN_POLL_SIZE`, `INITIAL_POLL_SIZE`)
3. Regenerate complete option docs
   - run `/Users/jracine/Development/crs/tools/nomad/generate_nomad_options_reference.sh`
   - confirm output at `/Users/jracine/Development/crs/inst/nomad/NOMAD_4_5_0_OPTIONS_REFERENCE.md` (or matching versioned filename)
   - verify `man/snomadr.Rd` still points users to the generated catalog
4. Benchmark and parity gates
   - run `run_nomad_smoke.R` and `run_nomad_bench.R`
   - compare against the latest accepted pre-upgrade baseline
   - include both fixed-seed and varying-seed comparisons
   - record mean/median elapsed deltas and objective/parameter drift
5. Exception-isolation gate (required when testing new solver families)
   - run `run_nomad_exception_isolation.R`
   - verify that an exception-only profile (for example `disco_opt`) does not break subsequent in-session runs of baseline/CS/NM profiles
   - keep isolated subprocess results as control
6. Important comparison caveat
   - `benchmarks/nomad/compare_prepost.R` merges by `case + seed_policy + seed`
   - repeated fixed seeds can cross-join replicates; for release decisions, use summary-level comparisons or add replicate index before merging
7. Optional heavy example gate (non-CRAN workflow)
   - use `/Users/jracine/Development/crs/man/runcrs` to enable running `\dontrun{}` examples
   - run tarball-first check from `/Users/jracine/Development`
   - restore with `/Users/jracine/Development/crs/man/dontruncrs` immediately after
   - remove disposable build artifacts (`src/*.o`, `src/*.so`) after check runs

## Exception-state guardrail

1. `solveNomadProblem()` in the embedded NOMAD C interface now performs exception-safe global cleanup on both success and failure paths.
2. This guardrail exists to prevent algorithm-flag/state contamination across sequential calls in one R session.
3. If contamination symptoms reappear in a future NOMAD update, run the exception-isolation harness first before trusting any solver-space benchmark conclusions.

## Package version/date workflow

`crs` release metadata is managed by:

- `/Users/jracine/Development/gen_crs` which runs:
  - `cd /Users/jracine/Development/crs`
  - `sh ./mkcrspkg.sh`

`mkcrspkg.sh` updates package metadata fields (for example `DESCRIPTION` date/version and `R/zzz.R` version string). These generated metadata updates are expected and safe to commit together with NOMAD integration changes.
