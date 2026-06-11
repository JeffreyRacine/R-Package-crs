# CRS Plot Route Inventory

Date: 2026-06-11

Scratch branch: `codex/crs-plot-modernization-scratch`

This inventory is the first tranche of the CRS plot modernization campaign. It
records the currently known public plot routes and the payload oracle that must
be preserved or intentionally superseded before a modern `np`-style renderer
owns the route.

## Routes

| Method | Current file | Current role | Initial oracle |
| --- | --- | --- | --- |
| `plot.crs` | `R/crs.R` | residual diagnostics by default; fitted slices with `mean=TRUE`; derivative slices with `deriv > 0`; rgl surface with `persp.rgl=TRUE` | `plot(..., mean=TRUE, plot.behavior="data")`; `plot(..., deriv=1, plot.behavior="data")`; existing rgl payload test |
| quantile `plot.crs` | `R/crs.R` | same branch as mean route when object has `tau` | payload generated from `plot(..., mean=TRUE, plot.behavior="data")` on a `tau` object |
| `plot.crsiv` | `R/crsiv.R` | univariate structural function / derivative route | sorted `xz[,1]`, `phi`, optional derivative matrix and interval matrices |
| `plot.crsivderiv` | `R/crsivderiv.R` | univariate derivative by default; structural function with `phi=TRUE` | sorted `xz[,1]`, `phi.prime` or `phi` |
| `plot.clsd` | `R/clsd.R` | density/distribution/derivative over training or evaluation grid | sorted `x`/`xer` paired with selected density/distribution vector |

## Initial Shadow Contract

The first internal payload builders live in `R/crs.plot.payload.R`. They are
not public API and they do not alter `plot()` defaults. They provide object-fed
payloads for future render engines:

- `.crs_plot_payload_regression()`
- `.crs_plot_payload_iv()`
- `.crs_plot_payload_iv_deriv()`
- `.crs_plot_payload_clsd()`

The regression payload builder has a `legacy=TRUE` mode solely for tests and
shadow comparisons against the existing `plot(..., plot.behavior="data")`
oracle.

## Guardrails

- No public CRS constructor or bandwidth search is called by the new payload
  builders.
- Public default plot behavior is unchanged in this tranche.
- Renderer behavior is unchanged in this tranche.
- Payload tests must compare statistical content, not only successful drawing.

## Scratch Notes

- During quantile plot-route proof, `crs(..., tau = 0.5, basis = "tensor")`
  with two continuous predictors failed during model fitting before any plot
  code was reached (`quantreg::rq()` received a non-scalar method). The same
  synthetic route fitted with `basis = "additive"`, and the opt-in plot surface
  matched legacy surface data. This observation is not treated as a plot-route
  regression gate in this tranche.
