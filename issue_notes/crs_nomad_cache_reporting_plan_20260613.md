# CRS NOMAD Cache Reporting Plan

Date: 2026-06-13

## Goal

Adopt the agreed user-facing NOMAD accounting style across CRS methods that expose NOMAD search summaries:

```text
Number of Function Evaluations: 252
NOMAD cache: 19,762 repeated point lookups avoided out of 20,014 (98.7%)
CRS restart cache: 14 repeated objective callbacks avoided out of 252 (5.6%)
```

The change is presentation-only. It must not change optimizer behavior, stored objective values, restart selection, cache semantics, or any numerical result.

## Current State

Most CRS NOMAD reporting is already centralized in `R/nomad.summary.R`:

- `summary.crs()` calls `.crs_nomad_summary_print()`.
- `print.crscv()`, `print.summary.crscv()`, and `summary.crscv()` call `.crs_nomad_summary_print()`.
- `summary.crsiv()`, `summary.crsivderiv()`, and `summary.clsd()` call `.crs_nomad_summary_print()`.
- `frscvNOMAD()` and `krscvNOMAD()` attach `nomad.summary` through `.crs_nomad_summary_from_solution()`.
- `clsd()` attaches `nomad.summary` for its NOMAD route.

The current CRS display already has the necessary counters:

- `num.feval`: blackbox/callback evaluations after NOMAD native cache.
- `cache.hits`: NOMAD native cache hits, cumulative over restarts.
- `total.evaluations`: NOMAD native point requests, cumulative over restarts.
- `callback.cache.hits`: CRS restart callback cache hits.
- `callback.cache.misses`: CRS restart callback cache misses.
- `callback.cache.size`: stored CRS restart callback cache entries.

## Implementation Plan

1. Update the central formatter in `R/nomad.summary.R`.
   - Change `.crs_nomad_native_cache_line()` to use the agreed `np`-style wording:
     `NOMAD cache: <hits> repeated point lookups avoided out of <requests> (<pct>%)`.
   - Keep the fallback for partial metadata conservative, but use wording that does not imply a denominator when none is available.

2. Update CRS restart-cache formatting in the same helper.
   - Change `.crs_nomad_restart_cache_line()` to print the compact agreed line when both hits and callback denominator are available:
     `CRS restart cache: <hits> repeated objective callbacks avoided out of <callbacks> (<pct>%)`.
   - Use `summary$callback.evaluations` as the denominator, because CRS restart-cache hits occur within the callback layer after NOMAD native caching.
   - Omit the restart-cache line when hits are zero, missing, or non-finite.
   - If hits are available but the denominator is missing, print a conservative fallback:
     `CRS restart cache: <hits> repeated objective callbacks avoided`.
   - Do not print `misses` or `stored points` in ordinary summaries; those counters remain internal diagnostics and can stay on the object.

3. Preserve object fields and internal accounting.
   - Do not rename `cache.hits`, `total.evaluations`, `callback.evaluations`, or callback-cache fields.
   - Do not change `.crs_nomad_restart_sweep()` accounting except if a denominator is missing from summaries despite already being available.
   - Do not alter `snomadr()` return contracts.

4. Sweep all NOMAD-facing summary paths.
   - Confirm the central helper covers:
     `summary.crs`, `print.crscv`, `print.summary.crscv`, `summary.crscv`,
     `summary.crsiv`, `summary.crsivderiv`, and `summary.clsd`.
   - Decide explicitly whether `print.snomadr()` should remain a low-level solver print method or adopt the same human-facing wording. Default: leave `print.snomadr()` unchanged unless it already exposes these high-level CRS summary lines.

5. Update tests.
   - Update `tests/testthat/test-nomad-summary.R` expected strings to the agreed wording.
   - Add a test for restart-cache percentage using `callback.evaluations` as denominator.
   - Add a zero-hit test to ensure no noisy cache line is printed.
   - Keep the live CRS summary smoke test to ensure the helper is still used in real model summaries.

6. Validate.
   - Run focused source-loaded tests:
     `tests/testthat/test-nomad-summary.R`.
   - Run adjacent source-loaded tests touching CRS search summaries:
     `tests/testthat/test-crs.R` and any NOMAD-native summary/counter tests.
   - Run full source-loaded `testthat` if focused tests are green.
   - Confirm with a small installed or source-loaded NOMAD example that output displays the agreed three-line style when both caches are nonzero.

## Risk Controls

- Keep this as a formatting-only tranche.
- Centralize changes in `R/nomad.summary.R`; avoid touching solver, restart, or CV objective code.
- Treat any change in `nomad.summary` object shape as out of scope.
- Preserve omitted-line behavior for missing or zero cache hits.
- Keep denominator checks strict: percentages are printed only when the denominator is finite, positive, and not smaller than the hit count.
- Do not duplicate the cache denominator logic in caller-specific summary methods.

## Acceptance Criteria

- All user-facing CRS NOMAD summary methods that use `.crs_nomad_summary_print()` display the agreed wording.
- NOMAD native cache percentage uses `cache.hits / total.evaluations`.
- CRS restart cache percentage uses `callback.cache.hits / callback.evaluations`.
- No cache line is printed for zero-hit caches.
- Numerical/model results are unchanged.
- Focused and adjacent tests are green.
