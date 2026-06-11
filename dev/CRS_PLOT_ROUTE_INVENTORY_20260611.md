# CRS Plot Route Inventory

Date: 2026-06-11

Scratch branch: `codex/crs-plot-modernization-scratch`

This inventory was the first tranche of the CRS plot modernization campaign.
It now also records the NP-interface direction that superseded the original
CRS-shaped opt-in route.

## Routes

| Method | Current file | Current role | Initial oracle |
| --- | --- | --- | --- |
| `plot.crs` | `R/crs.R` | fitted mean/quantile display by default; gradients via `gradients=TRUE`; rgl surface with `perspective = TRUE, renderer = "rgl"` | `plot(..., output = "data")`; `plot(..., gradients = TRUE, output = "data")`; existing rgl payload test |
| quantile `plot.crs` | `R/crs.R` | same branch as mean route when object has `tau` | payload generated from `plot(..., output ="data")` on a `tau` object |
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
shadow comparisons against the existing `plot(..., output = "data")`
oracle.

## Guardrails

- No public CRS constructor or bandwidth search is called by the new payload
  builders.
- Public `plot.crs()` behavior is intentionally changed from diagnostics to
  fitted-function display.
- Public plot spelling is moving to the modern NP interface.
- Payload tests must compare statistical content, not only successful drawing.

## Scratch Notes

- During quantile plot-route proof, `crs(..., tau = 0.5, basis = "tensor")`
  with two continuous predictors failed during model fitting before any plot
  code was reached (`quantreg::rq()` received a non-scalar method). The same
  synthetic route fitted with `basis = "additive"`, and the opt-in plot surface
  matched legacy surface data. This observation is not treated as a plot-route
  regression gate in this tranche.
