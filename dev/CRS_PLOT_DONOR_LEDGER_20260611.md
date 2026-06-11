# CRS Plot Donor-Code Ledger

Date: 2026-06-11

This ledger tracks plot code copied or adapted from `np` into the CRS scratch
plot modernization branch. It is deliberately conservative: only code that has
a CRS destination and a sentinel should be listed here.

| Donor | CRS destination | Mode | Sentinel |
| --- | --- | --- | --- |
| `np-master/R/np.plot.methods.R`: `np_boot_control()` | `R/crs.plot.helpers.R`: `.crs_boot_control()` | semantically adapted, CRS-internal name | `tests/testthat/test-plot-helpers.R` |
| `np-master/R/np.plot.methods.R`: `np_grid_control()` | `R/crs.plot.helpers.R`: `.crs_grid_control()` | semantically adapted, CRS-internal name | `tests/testthat/test-plot-helpers.R` |
| `np-master/R/np.plot.methods.R`: `np_render_control()` | `R/crs.plot.helpers.R`: `.crs_render_control()` | semantically adapted, CRS-internal name | `tests/testthat/test-plot-helpers.R` |
| `np-master/R/np.plot.methods.R`: public-dot validation helpers | `R/crs.plot.helpers.R`: `.crs_plot_validate_public_dots()` and support helpers | reduced CRS subset | `tests/testthat/test-plot-helpers.R` |
| `np-master/R/np.plot.methods.R`: `.np_plot_normalize_public_dots()` | `R/crs.plot.helpers.R`: `.crs_plot_normalize_public_dots()` | reduced CRS subset | `tests/testthat/test-plot-helpers.R` |
| `np-master/R/np.plot.helpers.R`: color/user-arg helpers | `R/crs.plot.helpers.R`: `.crs_plot_color()`, `.crs_plot_user_args()`, `.crs_plot_merge_user_args()` | reduced CRS subset | `tests/testthat/test-plot-helpers.R` |

Current status:

- copied helpers are internal-only;
- no public CRS plot default or renderer path uses them yet;
- no live CRS repository files are touched.

