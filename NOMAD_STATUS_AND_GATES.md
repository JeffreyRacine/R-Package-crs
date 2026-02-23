# NOMAD Status and Gates (`crs`)

Last updated: 2026-02-23

This file tracks dynamic findings, gate outcomes, and unresolved issues for the NOMAD4 migration.

## Scope completed so far

1. `crs` bridge migrated to embedded NOMAD `4.5.0`.
2. `src/Makevars` now builds against `src/nomad4_src`.
3. Legacy NOMAD3 tree removed from `src/nomad_src` (clean break).
4. NOMAD4 bridge behavior in `src/snomadr.cpp`:
   - direct NOMAD4 option names only
   - array/relative-option parsing
   - integer mesh/frame lower-bound enforcement
   - `MAX_EVAL` auto-derived from `MAX_BB_EVAL` when absent
   - `EPSILON` sanitization
5. Benchmark harness suite in place under `benchmarks/nomad/`.

## Primary migration lessons learned

1. NOMAD4 defaults are behaviorally different from NOMAD3 in ways that materially affect speed and selected parameters.
2. Without a `MAX_EVAL` cap, granular/integer runs can inflate evaluation budgets and produce severe runtime regressions.
3. A single global search-profile toggle is not safe:
   - settings that improve `frscvNOMAD`/`krscvNOMAD` can degrade `npglpreg` parity
   - settings that speed direct `snomadr` can alter solution/objective
4. Option handling must stay dimension-aware and type-aware (especially for integer/binary/categorical mapped inputs).

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

Clean-break rerun artifacts (after removing NOMAD3 tree and legacy aliases):

- `/tmp/crs_nomad_post_cleanbreak_20260223_1_raw.csv`
- `/tmp/crs_nomad_post_cleanbreak_20260223_1_summary.csv`
- `/tmp/crs_nomad_post_cleanbreak_20260223_1_parity.rds`
- `/tmp/crs_nomad_post_cleanbreak_20260223_1.log`
- `/tmp/crs_nomad_cleanbreak_vs_pre_summary_compare.csv`
- `/tmp/crs_nomad_cleanbreak_vs_pre_parity_compare.csv`
- `/tmp/crs_nomad_cleanbreak_vs_pre_minima_compare.csv`

## 2026-02-23 clean-break comparison vs NOMAD3 baseline

Reference pre-upgrade baseline:

- `/tmp/crs_nomad_pre_20260223_104449_raw.csv`
- `/tmp/crs_nomad_pre_20260223_104449_summary.csv`

Observed elapsed-time deltas (post vs pre, summary-level):

1. `frscvNOMAD`
   - fixed seed: mean `+665.85%`, median `+853.85%`
   - varying seed: mean `+954.84%`, median `+915.38%`
2. `krscvNOMAD`
   - fixed seed: mean `+578.64%`, median `+700.00%`
   - varying seed: mean `+981.82%`, median `+942.86%`
3. `npglpreg`
   - fixed seed: mean `-33.03%`, median `-34.36%`
   - varying seed: mean `+12.68%`, median `+33.23%`

Local-minima behavior across varying seeds (`npglpreg`, `frscvNOMAD`, `krscvNOMAD`):

1. `frscvNOMAD`: post objective better on `4/5` seeds
2. `krscvNOMAD`: post objective better on `3/5` seeds
3. `npglpreg`: post objective better on `5/5` seeds
   - best objective improved from `0.0421112191` (pre) to `0.036881` (post)
   - median objective improved from `0.0523705702` (pre) to `0.049461` (post)

Parity notes:

1. `frscvNOMAD` fixed-seed parameter signature changes: `0/5`
2. `krscvNOMAD` fixed-seed parameter signature changes: `0/5`
3. `npglpreg` fixed-seed parameter signature changes: `5/5` (drift remains)

Method note:

- `benchmarks/nomad/compare_prepost.R` can overcount fixed-seed repeats via cross-join.
- For this section, release-facing comparisons were taken from summary-level outputs and replicate-indexed parity/minima calculations.

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
