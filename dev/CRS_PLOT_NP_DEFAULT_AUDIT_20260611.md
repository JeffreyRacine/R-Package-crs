# CRS Plot NP Default Audit

Date: 2026-06-11

Purpose: track CRS plot defaults against the modern `np` plot layer so CRS
mean, quantile, IV, and CLSD plot methods converge on the same user-facing
behavior wherever the spline object model permits it.

## Donor Sources

- `np-master/R/np.plot.engine.rbandwidth.R`
- `np-master/R/np.plot.helpers.R`
- `np-master/R/np.plot.methods.R`
- CRS receiver files:
  - `R/crs.plot.helpers.R`
  - `R/crs.plot.engine.regression.R`
  - `R/crs.plot.engine.curves.R`
  - `R/crs.R`

## Shared Regression Defaults

| Control | NP default | CRS status |
|---|---:|---|
| `neval` | `50` | matched |
| `xq` | `0.5` | matched |
| scalar symmetric `xtrim` | `0.0` | matched for regression engines |
| `common.scale` | `TRUE` | matched |
| `perspective` | `TRUE` when exactly two continuous predictors | matched for supported CRS surfaces |
| `renderer` | `"base"` | matched |
| `view` | `"rotate"` | matched |
| `theta` | `0` | matched after default audit |
| `phi` | `20` for base; mapped to `-70` for rgl default view | matched after default audit |
| `plot.data.overlay` / `data_overlay` | `TRUE` | matched |
| `plot.rug` / `data_rug` | `FALSE` | matched |
| `plot.behavior` / `output` | `"plot"` | matched |
| `plot.par.mfrow` / `layout` | `TRUE` | matched for CRS regression slices |

## Shared Visual Roles

| Role | NP default | CRS status |
|---|---|---|
| fitted 1D line `type` | `"l"` | matched |
| fitted 1D line `col` | `par()$col` | matched after default audit |
| fitted 1D line `lwd` | `par()$lwd` | matched after default audit |
| fitted 1D line `lty` | `par()$lty` | matched after default audit |
| 1D plot `main` / `sub` | `""` | matched after default audit |
| overlay point `pch` | `20` | matched |
| overlay point `cex` | `0.5` | matched |
| overlay point color | viridis role `data_overlay`, alpha `0.35` | matched |
| support rug color | viridis role `support`, alpha `0.60` | matched |
| floor rug color | viridis role `support_floor`, alpha `0.55` | matched |
| floor/grid line widths | NP role widths | matched after default audit |
| surface palette | viridis facet colors, alpha `0.5` for base | matched |
| surface border | viridis role `surface_border`, role linewidth | matched after default audit |
| rgl surface alpha | `0.6` | matched |
| rgl overlay point size | `2` | matched |

## Interval Defaults

| Control | NP default | CRS status |
|---|---:|---|
| `errors` | `"none"` | matched |
| `band` | `"pmzsd"` / CRS internal `"standard"` | matched for supported routes |
| `alpha` | `0.05` | matched |
| bootstrap `B` | NP engine default `1999` | CRS public default remains `99` in current implementation; reassess separately because cost differs for spline refits |
| bootstrap method default | NP `"wild"` | CRS supports `"inid"` only and fails closed for others |
| bootstrap center | NP supports `"estimate"` / `"bias-corrected"` | CRS supports `"estimate"` only and fails closed otherwise |
| fitted 2D surface intervals | NP supports asymptotic/bootstrap | CRS supports asymptotic and inid-bootstrap after surface parity repair |
| derivative bootstrap intervals | NP supports where evaluator exists | CRS fails closed |

## Known Remaining Non-Identities

- CRS does not implement NP bootstrap methods `wild`, `fixed`, or `geom` for
  plot intervals; these fail closed rather than silently remapping to `inid`.
- CRS bootstrap default `B` is still `99`, not NP's engine default `1999`.
  Changing this would be a runtime/default-policy decision, not a visual
  parity cleanup.
- CRS surface gradients and bootstrap derivative intervals remain unsupported.
- `np_grid_control(xtrim=...)` endpoint-vector semantics should be reviewed in
  a separate focused tranche because the inherited NP regression engine itself
  still uses scalar symmetric trimming.
- `crsiv`, `crsivderiv`, and `clsd` use the CRS curve renderer and public
  NP-style argument grammar, but their statistical payloads are not identical
  to `npregiv` / `npregivderiv` objects. Further parity should be tested
  route-by-route with representative examples rather than by file copying.

## Acceptance Sentinels

- Exact user repro:
  `plot(crs(y ~ x + z, degree=c(5,1), segments=c(1,1), basis="auto", cv="none"), data_rug=TRUE, errors="asymptotic")`
- Focused plot tests:
  `plot-public-contract|plot-regression-shadow|plot-rgl|plot-helpers|plot-payload`
- Stress matrix:
  `dev/crs_plot_stress_matrix_20260611.R`

