# crs Pre-Release Smoke Remediation Plan

Date: 2026-04-29

## Post-Repair Status

Partly resolved in the live tree on 2026-04-29.

The clear demo-smoke defects were repaired without touching package algorithms:

- prompt-driven demos now use validated interactive input with noninteractive
  defaults and optional `CRS_DEMO_*` environment overrides;
- `radial_persp`, `radial_rgl`, and `sine_rgl` now define `nmulti` on both
  NOMAD and exhaustive branches;
- `nomad_factor` keeps its original large interactive default but uses a smaller
  noninteractive default.

Validation artifacts:

- `/Users/jracine/Development/tmp/pre_release_crs_demo_cleanup_20260429/logs/crsiv_default_after_patch.log`
- `/Users/jracine/Development/tmp/pre_release_crs_demo_cleanup_20260429/logs/radial_persp_default_after_patch.log`
- `/Users/jracine/Development/tmp/pre_release_crs_demo_cleanup_20260429/logs/radial_rgl_default_after_patch.log`
- `/Users/jracine/Development/tmp/pre_release_crs_demo_cleanup_20260429/logs/sine_rgl_default_after_patch.log`
- `/Users/jracine/Development/tmp/pre_release_crs_demo_cleanup_20260429/logs/nomad_factor_default_after_patch.log`
- `/Users/jracine/Development/tmp/pre_release_crs_demo_cleanup_20260429/logs/crs_check_as_cran_clean.log`

The tarball `R CMD check --as-cran --no-manual --no-vignettes
--ignore-vignettes` is green apart from the local `checkbashisms` warning and
the expected no-vignettes index note.

Deferred: `radial_constrained_test` still deserves a separate numerical-demo
diagnosis rather than a release-freeze band-aid.

## Scope

This plan covers issues found while exercising the submitted
`/Users/jracine/Development/crs_0.15-42.tar.gz` tarball through pre-release
smoke surfaces.

The primary release lanes are green:

- local `R CMD check --as-cran`: `1 WARNING` from missing local
  `checkbashisms`
- win-builder R-release: `OK`
- win-builder R-devel: `OK`
- full source test directory against private installed build: passed
- NOMAD benchmark smoke: passed

Primary artifacts:

- `/Users/jracine/Development/tmp/pre_release_smoke_20260429/test_logs/crs_source_full_testthat.log`
- `/Users/jracine/Development/tmp/pre_release_smoke_20260429/logs/crs_nomad_smoke.log`
- `/Users/jracine/Development/tmp/pre_release_smoke_20260429/results/crs_nomad_smoke_raw.csv`
- `/Users/jracine/Development/tmp/pre_release_smoke_20260429/results/crs_nomad_smoke_summary.csv`
- `/Users/jracine/Development/tmp/pre_release_smoke_20260429/results/crs_demo_status_full.csv`
- `/Users/jracine/Development/tmp/pre_release_smoke_20260429/demo_logs/crs`

## Issue CRS-1: Interactive Demos Fail In Noninteractive Smoke Runs

Affected demos:

- `crsiv`
- `crsiv_exog`
- `crsiv_exog_persp`
- `radial_persp`
- `radial_rgl`
- `sine_rgl`

Observed failure pattern:

- scripts call `readline()` for required choices;
- under noninteractive `Rscript`, empty input propagates to `NA`;
- subsequent `if (cv == "nomad")` checks fail with
  `missing value where TRUE/FALSE needed`.

## Issue CRS-2: `radial_constrained_test` Hits A Singular System

Observed failure:

```text
Error in solve.default(t(B) %*% B) :
  system is computationally singular
```

The same log also warns that the selected degree and segment counts equal the
search maxima, suggesting the demo configuration is pushing the example beyond
a stable smoke-demo regime.

## Issue CRS-3: `nomad_factor` Is Too Heavy For A Smoke Demo

Observed behavior:

- two successful `crs()` fits with `n = 10000` completed;
- the script exceeded the 300 second per-demo smoke timeout before finishing.

This is not evidence of estimator failure, but it is poor smoke-demo ergonomics.

## Highest-Standards Repair Strategy

1. Treat demos as user-facing examples, not tests to be gamed.
   They should run noninteractively with sensible defaults while preserving
   interactive prompts for users who launch them interactively.
2. For interactive demos:
   - introduce a small, shared demo input helper or a consistent local pattern;
   - use `interactive()` to decide whether to prompt;
   - in noninteractive mode, read optional environment variables first and then
     fall back to documented safe defaults;
   - validate user input before branching;
   - fail with a clear message only when an explicit invalid value is supplied.
3. For `radial_constrained_test`:
   - diagnose whether the singularity is due to data generation, selected basis
     dimension, constraint construction, or NOMAD selecting boundary complexity;
   - reduce the demo to a stable parameter regime if the purpose is
     illustration;
   - if the singularity exposes a real algorithmic weakness, repair the
     constrained-fit numerical contract at the appropriate helper and add a
     focused test.
4. For `nomad_factor`:
   - decide whether this belongs in `demo/` or in `benchmarks/nomad`;
   - if retained as a demo, reduce `n`, `nmulti`, or search budget so it is a
     reliable smoke example;
   - preserve the heavier version as a benchmark if it is scientifically useful.
5. Validation after candidate repair:
   - full source tests,
   - `R CMD check --as-cran` from tarball,
   - win-builder if native or compiled behavior changes,
   - noninteractive run of every `demo/*.R` with `R_DEFAULT_DEVICE=pdf` and
     `RGL_USE_NULL=TRUE`,
   - interactive spot-check for at least one prompt-preserving demo.

## Non-Issues From This Pass

- Bundled NOMAD compiler warnings remain accepted context from the current
  release hardening discussion.
- `cqrs`, `ggplot_cos`, `nomad_kernel`,
  `radial_constrained_first_partial`, `radial_constrained_mean`,
  `radial_constrained_second_partial`, and `spline` passed noninteractive demo
  smoke.

## Release Risk Classification

Low to medium.

The package check, win-builder, tests, and NOMAD smoke are green. The remaining
issues are demo quality and one unstable demo configuration. These should be
cleaned before release if demos are expected to serve as reliable user-facing
examples.
