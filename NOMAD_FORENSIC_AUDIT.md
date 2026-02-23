# NOMAD Forensic Audit (`crs`)

Audit date: 2026-02-23

Scope:

1. Active R integration layer:
   - `/Users/jracine/Development/crs/R`
2. Active native integration layer:
   - `/Users/jracine/Development/crs/src`
   - excluding upstream-vendored NOMAD source internals in `/Users/jracine/Development/crs/src/nomad4_src`

## Findings and actions

1. Legacy argument names in public R APIs (`min.poll.size.*`) remained in:
   - `crs.formula`
   - `frscvNOMAD`
   - `krscvNOMAD`
   - `npglpreg.formula`
   - `glpcvNOMAD`
2. Legacy-named default profile helper remained in `snomadr.R`:
   - `nomad4.compat.defaults()`
3. Stale/legacy comments remained in:
   - `/Users/jracine/Development/crs/R/crs.R`
   - `/Users/jracine/Development/crs/src/snomadr.cpp`
   - `/Users/jracine/Development/crs/src/crs_init.c`
4. Unused duplicate source tree remained:
   - `/Users/jracine/Development/crs/src/sgtelib_src`

## Remediation completed

1. Renamed API arguments from `min.poll.size.*` to `min.frame.size.*` across code and manuals:
   - `/Users/jracine/Development/crs/R/crs.R`
   - `/Users/jracine/Development/crs/R/frscvNOMAD.R`
   - `/Users/jracine/Development/crs/R/krscvNOMAD.R`
   - `/Users/jracine/Development/crs/R/np.regression.glp.R`
   - `/Users/jracine/Development/crs/man/crs.Rd`
   - `/Users/jracine/Development/crs/man/frscvNOMAD.Rd`
   - `/Users/jracine/Development/crs/man/krscvNOMAD.Rd`
   - `/Users/jracine/Development/crs/man/npglpreg.Rd`
2. Updated wrapper forwarding:
   - `crs.formula` now passes `min.frame.size.real` and `min.frame.size.integer` into `krscvNOMAD`.
3. Renamed NOMAD default helper to non-legacy naming:
   - `nomad4.compat.defaults()` -> `nomad4.mads.defaults()`
   - `merge.nomad4.compat.defaults()` -> `merge.nomad4.mads.defaults()`
4. Removed duplicate unused tree:
   - deleted `/Users/jracine/Development/crs/src/sgtelib_src`
5. Removed stale legacy/FIXME comments in active integration files.

## Validation

1. Build/install:
   - `/tmp/crs_install_after_forensic_legacy_cleanup_20260223.log`
2. API-formals checks:
   - `/tmp/crs_forensic_formals_check_20260223.out`
   - confirms `min.frame.size.*` present and `min.poll` absent for:
     - `crs.formula`
     - `krscvNOMAD`
     - `npglpreg.formula`
3. Runtime smoke:
   - `/tmp/crs_wrapper_minframe_smoke_20260223.out`
4. NOMAD smoke harness:
   - `/tmp/crs_nomad_smoke_after_forensic_legacy_cleanup_20260223_raw.csv`
   - `/tmp/crs_nomad_smoke_after_forensic_legacy_cleanup_20260223_summary.csv`
   - `/tmp/crs_nomad_smoke_after_forensic_legacy_cleanup_20260223_parity.rds`
5. Parity vs prior checkpoint:
   - `/tmp/crs_nomad_smoke_forensic_vs_wrapper_objdiff.csv`
   - result: `max_abs_objective_diff = 0` for all covered cases.

## Upstream caveat (vendor internals)

1. Vendored NOMAD4 upstream source includes compatibility internals (for example `ParametersNomad3.cpp`).
2. Attempting to exclude `ParametersNomad3.cpp` from build breaks link-time symbols (AllParameters vtable).
3. Current policy for `crs` integration:
   - keep upstream-required compilation units for NOMAD4 link correctness,
   - expose only NOMAD4 option names in `crs` APIs and docs,
   - keep legacy option names absent from active R and integration `src` code.
