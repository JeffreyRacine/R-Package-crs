# NOMAD Benchmark Harness (`crs`)

This folder provides strict pre/post parity and performance harnesses for NOMAD-related paths in `crs`.

Canonical result/status tracker:

- `/Users/jracine/Development/crs/NOMAD_STATUS_AND_GATES.md`

## Scripts

- `run_nomad_smoke.R`
  - small-`n`, fast validation run
  - writes:
    - `<prefix>_raw.csv`
    - `<prefix>_summary.csv`
    - `<prefix>_parity.rds`
- `run_nomad_bench.R`
  - larger run for more stable timing summaries
  - writes:
    - `<prefix>_raw.csv`
    - `<prefix>_summary.csv`
    - `<prefix>_parity.rds`
- `compare_prepost.R`
  - compares pre/post raw outputs and writes:
    - performance summary CSV (mean/median % elapsed change)
    - parity summary CSV (objective drift + parameter drift counts)
- `run_nomad_solver_space.R`
  - expanded solver/profile sweep for speed-vs-objective tradeoff mapping
  - writes:
    - `<prefix>_raw.csv`
    - `<prefix>_compare.csv`
    - `<prefix>_agg.csv`
    - `<prefix>_strict_safe.csv`
- `run_nomad_solver_space_extended.R`
  - deeper sweep adding MADS search-method toggles (`SIMPLE_LINE_SEARCH`, `SPECULATIVE_SEARCH`, `NM_SEARCH`, `VNS_MADS_SEARCH`) plus standalone solver modes
  - writes:
    - `<prefix>_raw.csv`
    - `<prefix>_compare.csv`
    - `<prefix>_agg.csv`
    - `<prefix>_strict_safe.csv`
- `run_nomad_mads_deep_space.R`
  - MADS-only deep sweep for direction/search/restart/budget controls
  - writes:
    - `<prefix>_raw.csv`
    - `<prefix>_compare.csv`
    - `<prefix>_agg.csv`
    - `<prefix>_strict_safe.csv`
    - `<prefix>_relaxed_safe.csv`
- `run_nomad_kr_timebudget.R`
  - `krscvNOMAD`-focused speed/accuracy frontier sweep including time-budget reinvestment profiles
  - writes:
    - `<prefix>_raw.csv`
    - `<prefix>_compare.csv`
    - `<prefix>_agg.csv`
    - `<prefix>_strict_safe.csv`
- `run_nomad_dimmix_tradeoff.R`
  - mixed-dimension tradeoff sweep for `(p_cont, p_cat) = (1,1), (2,2), (3,3)` at fixed `n`
  - compares baseline vs aggressive profiles for:
    - `frscvNOMAD`
    - `krscvNOMAD`
    - `npglpreg`
  - writes:
    - `<prefix>_raw.csv`
    - `<prefix>_compare.csv`
    - `<prefix>_agg_by_dim.csv`
    - `<prefix>_agg_overall.csv`
- `run_nomad_exception_isolation.R`
  - verifies that exception-only profiles do not poison later calls in the same R session
  - compares in-session sequence vs isolated subprocess calls
  - writes:
    - `<prefix>_raw.csv`
    - `<prefix>_summary.csv`

## Coverage

Each suite runs:

- `frscvNOMAD`
- `krscvNOMAD`
- `npglpreg` (NOMAD-driven CV path)
- `snomadr_basic_lib` (direct `snomadr` case extracted from `man/snomadr.Rd`)

with both:

- fixed-seed repetitions
- varying-seed repetitions

## Example commands

```bash
cd /Users/jracine/Development/crs

Rscript benchmarks/nomad/run_nomad_smoke.R /tmp/crs_nomad_pre_smoke
Rscript benchmarks/nomad/run_nomad_bench.R /tmp/crs_nomad_pre_bench

# ...after changes...
Rscript benchmarks/nomad/run_nomad_smoke.R /tmp/crs_nomad_post_smoke
Rscript benchmarks/nomad/run_nomad_bench.R /tmp/crs_nomad_post_bench

Rscript benchmarks/nomad/compare_prepost.R \
  /tmp/crs_nomad_pre_bench_raw.csv \
  /tmp/crs_nomad_post_bench_raw.csv \
  /tmp/crs_nomad_prepost_perf_summary.csv \
  /tmp/crs_nomad_prepost_parity_summary.csv

# expanded solver-space sweep (second arg is n; optional)
Rscript benchmarks/nomad/run_nomad_solver_space.R \
  /tmp/crs_nomad_solver_space \
  180

# deeper solver-space sweep (second arg is n; optional)
Rscript benchmarks/nomad/run_nomad_solver_space_extended.R \
  /tmp/crs_nomad_solver_space_ext \
  180

# MADS-only deep sweep (second arg is n; optional)
Rscript benchmarks/nomad/run_nomad_mads_deep_space.R \
  /tmp/crs_nomad_mads_deep_space \
  180

# kr time-budget frontier sweep (second arg is n; optional)
Rscript benchmarks/nomad/run_nomad_kr_timebudget.R \
  /tmp/crs_nomad_kr_timebudget \
  300

# mixed-dimension tradeoff sweep (second arg is n; optional)
Rscript benchmarks/nomad/run_nomad_dimmix_tradeoff.R \
  /tmp/crs_nomad_dimmix_tradeoff \
  100

# exception-isolation gate
Rscript benchmarks/nomad/run_nomad_exception_isolation.R \
  /tmp/crs_nomad_exception_isolation \
  180
```
