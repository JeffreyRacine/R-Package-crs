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
```
