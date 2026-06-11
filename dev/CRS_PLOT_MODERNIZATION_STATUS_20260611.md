# CRS Plot Modernization Scratch Status

Date: 2026-06-11

Branch: `codex/crs-plot-modernization-scratch`

Live repo mutation: none.

## Candidate State

This scratch branch contains an opt-in modern plot candidate. It does not
change legacy defaults.

Opt-in selector:

- `plot.view = "fit"`

Legacy defaults retained:

- `plot.crs()` still enters the diagnostic path.
- `plot.crs(mean = TRUE)` still enters the legacy fitted/mean path.
- `plot.crsiv()`, `plot.crsivderiv()`, and `plot.clsd()` still enter their
  legacy paths unless `plot.view = "fit"` is supplied.

## Implemented Routes

- Object-fed regression/quantile payloads for `crs` objects.
- Object-fed IV payloads for `crsiv`.
- Object-fed IV derivative-object payloads for `crsivderiv`.
- Object-fed CLSD density/distribution/derivative payloads.
- Minimal CRS-owned plot donor helpers adapted from `np`.
- Opt-in 1D CRS fitted/quantile route with asymptotic interval parity.
- Opt-in 2D CRS point-estimate surface route for exactly two continuous
  predictors.
- Opt-in `crsiv`, `crsivderiv`, and `clsd` curve routes with `output="data"`.

## Fail-Closed Boundaries

- Bootstrap intervals are not enabled through the opt-in CRS route yet.
- 2D interval surfaces are not enabled yet.
- Derivative CRS plots are not enabled through the modern regression payload
  route yet.
- Contradictory `persp.rgl=TRUE` and `renderer!="rgl"` errors.
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

Package check:

- Command:
  `R CMD check --no-manual --no-build-vignettes --library=/Users/jracine/Development/tmp/crs_plot_modernization_20260611/Rlib /Users/jracine/Development/tmp/crs_plot_modernization_20260611/check/crs_0.15-45.tar.gz`
- Result: completed.
- Warnings: two vignette-output warnings caused by checking with
  `--no-build-vignettes` and no `inst/doc`; examples and vignette R code
  passed.

## Scratch Observation

A synthetic two-continuous `crs(..., tau = 0.5, basis = "tensor")` fit failed
during model construction before plot code was reached. The additive analogue
fit successfully and matched legacy surface data through the opt-in plot route.
This is recorded as a non-plot observation, not as a plot-route regression.

## Remaining Promotion Decision

The next high-risk decision is whether and when to switch the public
`plot.crs()` default from diagnostics to fitted function display. That should
be a separate tranche with documentation and NEWS updates because it changes
user-facing behavior even though the opt-in route is now proven in scratch.
