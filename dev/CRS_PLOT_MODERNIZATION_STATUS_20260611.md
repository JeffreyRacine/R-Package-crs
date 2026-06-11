# CRS Plot Modernization Scratch Status

Date: 2026-06-11

Branch: `codex/crs-plot-modernization-scratch`

Live repo mutation: none.

## Candidate State

This scratch branch contains an NP-interface plot candidate for CRS objects.
`plot.crs()` now displays fitted function surfaces by default and uses
NP-style public plot arguments.

Selector:

- `errors`
- `gradients`
- `gradient_order`
- `output`
- `perspective`
- `renderer`
- `neval`

Default behavior:

- `plot.crs()` enters the fitted mean/quantile path.
- The old lm-style diagnostic panel has been removed from `plot.crs()`.
- `plot.crsiv()`, `plot.crsivderiv()`, and `plot.clsd()` now enter their
  modern curve routes directly.

## Implemented Routes

- Object-fed regression/quantile payloads for `crs` objects.
- Object-fed IV payloads for `crsiv`.
- Object-fed IV derivative-object payloads for `crsivderiv`.
- Object-fed CLSD density/distribution/derivative payloads.
- Minimal CRS-owned plot donor helpers adapted from `np`.
- Default 1D CRS fitted/quantile route with asymptotic interval parity and
  bootstrap mean intervals.
- Default CRS derivative/gradient data and render route with asymptotic
  intervals.
- Default 2D CRS point-estimate surface route for exactly two continuous
  predictors.
- Direct `crsiv`, `crsivderiv`, and `clsd` curve routes with `output="data"`.

## Fail-Closed Boundaries

- Bootstrap derivative intervals are not enabled and fail closed.
- 2D interval surfaces are not enabled yet.
- Contradictory `perspective = TRUE, renderer = "rgl"` and `renderer!="rgl"` errors.
- Unsupported surface predictor shapes error through the payload builder.

## Checkpoints

- `ccca45c` - CRS plot payload shadow scaffold.
- `2dcc354` - minimal CRS plot donor helpers.
- `23b6ee1` - shadow 1D regression route.
- `44ae1a0` - opt-in CRS fit plot route.
- `d8f607a` - opt-in CRS surface fit route.
- `2e0654c` - quantile fit plot parity proof.
- `5001062` - fit interval plot parity proof.
- `cc8e61b` - opt-in modern curve plot routes.
- `83962cb` - switch `plot.crs()` default to fitted displays and remove the
  diagnostic monolith.
- `d6a4613` - adopt NP-style CRS plot interface in scratch.

## Installed Proof

Passed from the scratch-installed package:

- `test-plot-curves-modern.R`
- `test-plot-regression-shadow.R`
- `test-plot-payload.R`
- `test-plot-helpers.R`
- `test-plot-rgl.R`
- `test-crsiv.R`
- `test-crsivderiv.R`
- `test-clsd.R`
- `test-crs.R`

Additional default-switch proof:

- `R CMD check --no-manual --no-build-vignettes` completed on the scratch
  tarball after the default switch.

Package check:

- Command:
  `R CMD check --no-manual --library=/Users/jracine/Development/tmp/crs_plot_modernization_20260611/Rlib /Users/jracine/Development/tmp/crs_plot_modernization_20260611/check_stress/crs_0.15-45.tar.gz`
- Result: `Status: OK`.

## Stress Gate

The archived `np` plot closeout matrices were adapted into a CRS scratch stress
runner:

- Script:
  `dev/crs_plot_stress_matrix_20260611.R`
- Result directory:
  `/Users/jracine/Development/tmp/crs_plot_modernization_20260611/plot_stress_20260611`
- Result: 30/30 cases passed with `B = 5` and `neval = 7`.

Coverage includes:

- CRS fitted mean data/render routes.
- CRS asymptotic and inid-bootstrap interval payloads.
- CRS gradient payloads and asymptotic gradient intervals.
- CRS 2D fitted surfaces.
- CRS quantile fitted and gradient routes.
- `crsiv`, `crsivderiv`, and `clsd` data/render routes.
- Fail-fast cases for legacy CRS-only plot arguments, unsupported interval
  combinations, unsupported bootstrap combinations, and unknown arguments.

The first stress run exposed a real public-adapter bug: `plot.crs()` evaluated
unknown `...` promises before rejecting them. The adapter now validates the
unevaluated call object before forcing `list(...)`, and
`test-plot-regression-shadow.R` contains a focused regression test for this
contract.

## Scratch Observation

A synthetic two-continuous `crs(..., tau = 0.5, basis = "tensor")` fit failed
during model construction before plot code was reached. The additive analogue
fit successfully and matched legacy surface data through the opt-in plot route.
This is recorded as a non-plot observation, not as a plot-route regression.

## Remaining Promotion Work

- Live repo promotion still requires explicit approval.
- NEWS/CHANGELOG language should call out the intentional `plot.crs()`
  default change from diagnostics to fitted function display.
- 2D interval surfaces and bootstrap derivative intervals remain future
  fail-closed extension points.
