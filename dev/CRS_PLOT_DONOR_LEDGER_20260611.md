# CRS Plot Donor-Code Ledger

Date: 2026-06-11

This ledger tracks plot code copied or adapted from `np` into the CRS scratch
plot modernization branch. It is deliberately conservative: only code that has
a CRS destination and a sentinel should be listed here.

| Donor | CRS destination | Mode | Sentinel |
| --- | --- | --- | --- |
| `np-master/R/np.plot.methods.R`: `np_boot_control()` | `R/crs.plot.helpers.R`: `np_boot_control()` and `crs_boot_control()` alias | copied/adapted with NP class accepted by CRS plots | `tests/testthat/test-plot-helpers.R` |
| `np-master/R/np.plot.methods.R`: `np_grid_control()` | `R/crs.plot.helpers.R`: `np_grid_control()` and `crs_grid_control()` alias | copied/adapted with NP class accepted by CRS plots | `tests/testthat/test-plot-helpers.R` |
| `np-master/R/np.plot.methods.R`: `np_render_control()` | `R/crs.plot.helpers.R`: `np_render_control()` and `crs_render_control()` alias | copied/adapted with NP class accepted by CRS plots | `tests/testthat/test-plot-helpers.R` |
| `np-master/R/np.plot.methods.R`: public-dot validation helpers | `R/crs.plot.helpers.R`: `.crs_plot_validate_public_dots()` and support helpers | NP-shaped CRS adapter | `tests/testthat/test-plot-helpers.R` |
| `np-master/R/np.plot.methods.R`: `.np_plot_normalize_public_dots()` | `R/crs.plot.helpers.R`: `.crs_plot_normalize_public_dots()` | NP-shaped CRS adapter | `tests/testthat/test-plot-helpers.R` |
| `np-master/R/np.plot.helpers.R`: color/user-arg helpers | `R/crs.plot.helpers.R`: `.crs_plot_color()`, `.crs_plot_user_args()`, `.crs_plot_merge_user_args()` | reduced CRS subset | `tests/testthat/test-plot-helpers.R` |

Current status:

- NP-named control helpers are exported by the scratch package;
- public `plot.crs()` routes through the CRS normalizer and accepts NP-shaped
  controls while using `predict.crs()`/spline helpers as the evaluation backend;
- no live CRS repository files are touched.
