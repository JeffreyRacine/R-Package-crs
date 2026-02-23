# NOMAD Status and Gates (`crs`)

Last updated: 2026-02-23

This file tracks dynamic findings, gate outcomes, and unresolved issues for the NOMAD4 migration.

## Scope completed so far

1. `crs` bridge migrated to embedded NOMAD `4.5.0`.
2. `src/Makevars` now builds against `src/nomad4_src`.
3. Legacy option compatibility shims added in `src/snomadr.cpp`:
   - `MIN_POLL_SIZE` -> `MIN_FRAME_SIZE`
   - `INITIAL_POLL_SIZE` -> `INITIAL_FRAME_SIZE`
   - array/relative-option parsing
   - integer mesh/frame lower-bound enforcement
   - `MAX_EVAL` auto-derived from `MAX_BB_EVAL` when absent
   - `EPSILON` sanitization
4. Benchmark harness suite in place under `benchmarks/nomad/`.

## Primary migration lessons learned

1. NOMAD4 defaults are behaviorally different from NOMAD3 in ways that materially affect speed and selected parameters.
2. Without a `MAX_EVAL` cap, granular/integer runs can inflate evaluation budgets and produce severe runtime regressions.
3. A single global search-profile toggle is not safe:
   - settings that improve `frscvNOMAD`/`krscvNOMAD` can degrade `npglpreg` parity
   - settings that speed direct `snomadr` can alter solution/objective
4. Option translation must stay dimension-aware and type-aware (especially for integer/binary/categorical mapped inputs).

## Latest gate artifacts

Baseline (pre-upgrade reference):

- `/tmp/crs_nomad_pre_20260223_104449_raw.csv`
- `/tmp/crs_nomad_pre_20260223_104449_summary.csv`

Latest post candidate:

- `/tmp/crs_nomad_post_20260223_104449_v2_raw.csv`
- `/tmp/crs_nomad_post_20260223_104449_v2_summary.csv`

Latest pre/post comparisons:

- `/tmp/crs_nomad_perf_20260223_104449_v2.csv`
- `/tmp/crs_nomad_parity_20260223_104449_v2.csv`

Supporting logs:

- `/tmp/crs_nomad_gate_20260223_104449.log`
- `/tmp/crs_nomad_gate_20260223_104449_v2.log`

Recent smoke after clean rebuild:

- `/tmp/crs_nomad_post_smoke_final_raw.csv`
- `/tmp/crs_nomad_post_smoke_final_summary.csv`
- `/tmp/crs_nomad_post_smoke_final_parity.rds`
- `/tmp/crs_nomad_smoke_final.log`

## Current gate state

1. Build/install: passing.
2. Runtime compatibility: passing for core paths (`frscvNOMAD`, `krscvNOMAD`, `npglpreg`, `snomadr`).
3. Strict parity/performance gates: not yet fully passing.

Outstanding gaps:

1. Fixed-seed performance regressions remain in `frscvNOMAD` and `krscvNOMAD` vs NOMAD3 baseline.
2. `npglpreg` still shows parameter drift under strict comparisons.
3. Direct `snomadr` fixed-seed behavior is sensitive to search-profile changes.

## Practical forward plan

1. Keep bridge contract stable and continue tuning at option-profile level.
2. Prefer path-aware tuning over global search-disable flags.
3. Re-run strict pre/post gates after each tuning checkpoint and append artifacts here.
4. Treat this file as the only canonical place for migration status and gate outcomes.
