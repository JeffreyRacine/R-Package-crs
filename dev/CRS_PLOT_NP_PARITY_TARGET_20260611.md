# CRS Plot NP Parity Target

This scratch tranche supersedes the earlier CRS-shaped plot modernization
candidate. The promotion target is wholesale adoption of the modern `np` plot
interface and behavior, adapted only where CRS uses spline evaluators rather
than local-kernel bandwidth objects.

## Non-Negotiable Target

- `plot.crs()` mirrors `plot.npregression()` / `plot.qregression()`:
  fitted-function display by default, `gradients=TRUE` for derivative displays,
  `errors=` for intervals, `output=` for return mode, `perspective=`,
  `renderer=`, `view=`, `neval=`, `grid_control=`, `boot_control=`, and
  `render_control=`.
- `plot.crs(..., tau=)` / fitted quantile objects mirror `plot.npqreg()`
  semantics where the CRS object is a conditional-quantile fit.
- `plot.crsiv()` mirrors `plot.npregiv()` to the extent the current CRS IV
  object contains the same structural-function and derivative payloads.
- `plot.crsivderiv()` mirrors `plot.npregivderiv()`.
- `plot.clsd()` should use the same modern output/render controls rather than
  legacy CRS-only plot switches.

## Implementation Architecture

1. Treat `predict.crs()` as the CRS analogue of the NP plot hat/evaluation
   helpers. It already accepts `newdata`, supports `deriv`, and returns fitted
   values plus interval and derivative attributes.
2. Copy/adapt NP's public argument grammar and validation helpers into CRS:
   `errors`, `band`, `alpha`, `bootstrap`, `B`, `center`, `output`,
   `data_overlay`, `data_rug`, `layout`, `legend`, `gradient`, `gradients`,
   `gradient_order`, `common_scale`, `renderer`, `neval`, `perspective`,
   `view`, `np_boot_control`, `np_grid_control`, and `np_render_control`.
3. Keep CRS standalone. Do not call `np::` or `np:::` from package internals.
4. Fail closed for NP options until the CRS evaluator genuinely implements
   them. No silent downgrades.
5. Port NP output-object shape where user code depends on it; data-frame-only
   payloads are an intermediate scratch state, not final parity.

## Current Parity Gaps

- CRS `output="data"` currently returns simple slice/surface data frames in the
  scratch route; NP generally returns estimator-like payload objects. This must
  be reconciled before promotion.
- CRS bootstrap intervals currently implement only the existing resampling
  route; NP bootstrap choices such as `wild`, `fixed`, and `geom` must either
  be implemented or rejected explicitly by route.
- CRS fitted surface intervals now work for asymptotic and inid-bootstrap
  mean/quantile surfaces, including base and `rgl` renderers. Surface gradients
  remain unsupported and must fail closed until a CRS evaluator route exists.
- CRS IV and CLSD routes still require full NP-argument cleanup after the main
  `plot.crs()` adapter is stabilized.
- `np_grid_control(xtrim=...)` parity should be revisited separately: the
  modern constructor follows endpoint-quantile validation, while the inherited
  NP regression plot engines still use scalar symmetric trimming internally.

## Acceptance Gates

- Installed package tests must demonstrate NP spelling for all public examples:
  `errors`, `gradients`, `gradient_order`, `output`,
  `perspective`, `renderer`, `neval`, `view`, and public CRS control helpers.
- Legacy CRS examples such as `plot(model, ci=TRUE)` and
  `plot(model, deriv=1)` must not appear in promoted docs/tests.
- Public legacy arguments must either be absent or fail with explicit
  remediation pointing to the NP argument.
- Plot parity sentinels must cover mean, quantile, gradients, 2D surface,
  bootstrap intervals, IV structural function, IV derivative, and CLSD.
- Surface sentinels must include `data_rug`, `data_overlay`, base `persp`,
  `renderer="rgl"`, asymptotic intervals, and bootstrap intervals.
