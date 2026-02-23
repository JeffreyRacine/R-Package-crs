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

## Current strategic baseline

1. `crs` is wired to embedded NOMAD `4.5.0`.
2. `crs` keeps a NOMAD4-only source layout (legacy NOMAD3 tree removed).
3. Near-term goal remains: stabilize NOMAD4 under strict parity/performance gates for existing `crs` behavior.

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
5. Important comparison caveat
   - `benchmarks/nomad/compare_prepost.R` merges by `case + seed_policy + seed`
   - repeated fixed seeds can cross-join replicates; for release decisions, use summary-level comparisons or add replicate index before merging
6. Optional heavy example gate (non-CRAN workflow)
   - use `/Users/jracine/Development/crs/man/runcrs` to enable running `\dontrun{}` examples
   - run tarball-first check from `/Users/jracine/Development`
   - restore with `/Users/jracine/Development/crs/man/dontruncrs` immediately after
   - remove disposable build artifacts (`src/*.o`, `src/*.so`) after check runs
