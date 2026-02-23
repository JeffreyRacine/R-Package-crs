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

## Current strategic baseline

1. `crs` is wired to embedded NOMAD `4.5.0`.
2. Legacy NOMAD3 source is retained in-tree for reference only.
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
