# NOMAD Status and Gates (`crs`)

Last updated: 2026-02-23

This file tracks dynamic findings, gate outcomes, and unresolved issues for the NOMAD4 migration.

## Scope completed so far

1. `crs` bridge migrated to embedded NOMAD `4.5.0`.
2. `src/Makevars` now builds against `src/nomad4_src`.
3. Legacy NOMAD3 tree removed from `src/nomad_src` (clean break).
4. NOMAD4 bridge behavior in `src/snomadr.cpp`:
   - direct NOMAD4 option names only
   - array/relative-option parsing
   - integer mesh/frame lower-bound enforcement
   - `MAX_EVAL` auto-derived from `MAX_BB_EVAL` when absent
   - `EPSILON` sanitization
5. Benchmark harness suite in place under `benchmarks/nomad/`.

## Primary migration lessons learned

1. NOMAD4 defaults are behaviorally different from NOMAD3 in ways that materially affect speed and selected parameters.
2. Without a `MAX_EVAL` cap, granular/integer runs can inflate evaluation budgets and produce severe runtime regressions.
3. A single global search-profile toggle is not safe:
   - settings that improve `frscvNOMAD`/`krscvNOMAD` can degrade `npglpreg` parity
   - settings that speed direct `snomadr` can alter solution/objective
4. Option handling must stay dimension-aware and type-aware (especially for integer/binary/categorical mapped inputs).

## Latest gate artifacts

Baseline (pre-upgrade reference):

- `/tmp/crs_nomad_pre_20260223_104449_raw.csv`
- `/tmp/crs_nomad_pre_20260223_104449_summary.csv`

Latest post candidate:

- `/tmp/crs_nomad_post_20260223_104449_v2_raw.csv`
- `/tmp/crs_nomad_post_20260223_104449_v2_summary.csv`

Latest pre/post comparisons:

- `/tmp/crs_nomad_perf_20260223_104449_v2.csv`
- `/tmp/crs_nomad_parity_20260223_104449_v2.csv`

Supporting logs:

- `/tmp/crs_nomad_gate_20260223_104449.log`
- `/tmp/crs_nomad_gate_20260223_104449_v2.log`

Recent smoke after clean rebuild:

- `/tmp/crs_nomad_post_smoke_final_raw.csv`
- `/tmp/crs_nomad_post_smoke_final_summary.csv`
- `/tmp/crs_nomad_post_smoke_final_parity.rds`
- `/tmp/crs_nomad_smoke_final.log`

Clean-break rerun artifacts (after removing NOMAD3 tree and legacy aliases):

- `/tmp/crs_nomad_post_cleanbreak_20260223_1_raw.csv`
- `/tmp/crs_nomad_post_cleanbreak_20260223_1_summary.csv`
- `/tmp/crs_nomad_post_cleanbreak_20260223_1_parity.rds`
- `/tmp/crs_nomad_post_cleanbreak_20260223_1.log`
- `/tmp/crs_nomad_cleanbreak_vs_pre_summary_compare.csv`
- `/tmp/crs_nomad_cleanbreak_vs_pre_parity_compare.csv`
- `/tmp/crs_nomad_cleanbreak_vs_pre_minima_compare.csv`

Compatibility-profile rerun artifacts (explicit NOMAD3.9.1-like defaults via NOMAD4 names):

- `/tmp/crs_install_compat_profile_20260223.log`
- `/tmp/crs_nomad_post_compatdefaults_20260223_1_raw.csv`
- `/tmp/crs_nomad_post_compatdefaults_20260223_1_summary.csv`
- `/tmp/crs_nomad_post_compatdefaults_20260223_1_parity.rds`
- `/tmp/crs_nomad_post_compatdefaults_20260223_1.log`
- `/tmp/crs_nomad_compat_vs_pre_20260223_1_summary_compare.csv`
- `/tmp/crs_nomad_compat_vs_pre_20260223_1_parity_compare.csv`
- `/tmp/crs_nomad_compat_vs_pre_20260223_1_minima_compare.csv`
- `/tmp/crs_nomad_compat_vs_prevpost_20260223_1_summary_compare.csv`

Compatibility-profile refresh rerun (same profile, fresh install/build):

- `/tmp/crs_install_compatprofile_refresh_20260223.log`
- `/tmp/crs_nomad_post_compatdefaults_20260223_2_raw.csv`
- `/tmp/crs_nomad_post_compatdefaults_20260223_2_summary.csv`
- `/tmp/crs_nomad_post_compatdefaults_20260223_2_parity.rds`
- `/tmp/crs_nomad_post_compatdefaults_20260223_2.log`
- `/tmp/crs_nomad_compat_vs_pre_20260223_2_summary_compare.csv`
- `/tmp/crs_nomad_compat_vs_pre_20260223_2_parity_compare.csv`
- `/tmp/crs_nomad_compat_vs_pre_20260223_2_minima_compare.csv`
- `/tmp/crs_nomad_compat_vs_prevpost_20260223_2_summary_compare.csv`
- `/tmp/crs_nomad_compat_vs_prevpost_20260223_2_parity_compare.csv`
- `/tmp/crs_nomad_compat_vs_prevpost_20260223_2_minima_compare.csv`

Path-specific profile experiments (MADS search-option tuning):

- `/tmp/crs_nomad_profile_scan_20260223_raw.csv`
- `/tmp/crs_nomad_profile_scan_20260223_summary.csv`
- `/tmp/crs_nomad_profile_scan_20260223_vs_base.csv`
- `/tmp/crs_nomad_path_tuned_confirm_20260223_raw.csv`
- `/tmp/crs_nomad_path_tuned_confirm_20260223_summary.csv`
- `/tmp/crs_nomad_path_tuned_confirm_20260223_compare.csv`
- `/tmp/crs_nomad_path_profiles_compare_20260223_raw.csv`
- `/tmp/crs_nomad_path_profiles_compare_20260223_summary.csv`
- `/tmp/crs_nomad_path_profiles_compare_20260223_vs_base.csv`
- `/tmp/crs_nomad_path_no_quad_vs_pre_20260223_summary.csv`
- `/tmp/crs_nomad_path_no_quad_vs_pre_20260223_minima.csv`
- `/tmp/crs_nomad_path_cs_vs_pre_20260223_summary.csv`
- `/tmp/crs_nomad_path_cs_vs_pre_20260223_minima.csv`

Standalone solver-mode scan (NOMAD4 algorithm toggles):

- `/tmp/crs_nomad_solver_modes_scan_20260223_raw.csv`
- `/tmp/crs_nomad_solver_modes_scan_20260223_summary.csv`
- `/tmp/crs_nomad_solver_modes_scan_20260223.log`

Post-implementation benchmark (path defaults active in package code):

- `/tmp/crs_install_pathdefaults_20260223.log`
- `/tmp/crs_nomad_post_pathdefaults_20260223_1_raw.csv`
- `/tmp/crs_nomad_post_pathdefaults_20260223_1_summary.csv`
- `/tmp/crs_nomad_post_pathdefaults_20260223_1_parity.rds`
- `/tmp/crs_nomad_post_pathdefaults_20260223_1.log`
- `/tmp/crs_nomad_pathdefaults_vs_pre_20260223_1_summary_compare.csv`
- `/tmp/crs_nomad_pathdefaults_vs_pre_20260223_1_parity_compare.csv`
- `/tmp/crs_nomad_pathdefaults_vs_pre_20260223_1_minima_compare.csv`
- `/tmp/crs_nomad_pathdefaults_vs_compat_20260223_1_summary_compare.csv`
- `/tmp/crs_nomad_pathdefaults_vs_compat_20260223_1_parity_compare.csv`
- `/tmp/crs_nomad_pathdefaults_vs_compat_20260223_1_minima_compare.csv`

Extended solver-space follow-up (search-method toggles plus solver modes):

- `/tmp/crs_nomad_solver_space_ext_20260223_raw.csv`
- `/tmp/crs_nomad_solver_space_ext_20260223_compare.csv`
- `/tmp/crs_nomad_solver_space_ext_20260223_agg.csv`
- `/tmp/crs_nomad_solver_space_ext_20260223_strict_safe.csv`
- `/tmp/crs_nomad_solver_space_ext_20260223.log`

MADS-only deep sweeps (direction/search/restart/budget controls):

- `/tmp/crs_nomad_mads_deep_space_20260223_raw.csv`
- `/tmp/crs_nomad_mads_deep_space_20260223_compare.csv`
- `/tmp/crs_nomad_mads_deep_space_20260223_agg.csv`
- `/tmp/crs_nomad_mads_deep_space_20260223_strict_safe.csv`
- `/tmp/crs_nomad_mads_deep_space_20260223_relaxed_safe.csv`
- `/tmp/crs_nomad_mads_deep_space_20260223.log`
- `/tmp/crs_nomad_mads_deep_space_20260223_n300_raw.csv`
- `/tmp/crs_nomad_mads_deep_space_20260223_n300_compare.csv`
- `/tmp/crs_nomad_mads_deep_space_20260223_n300_agg.csv`
- `/tmp/crs_nomad_mads_deep_space_20260223_n300_strict_safe.csv`
- `/tmp/crs_nomad_mads_deep_space_20260223_n300_relaxed_safe.csv`
- `/tmp/crs_nomad_mads_deep_space_20260223_n300.log`

Post-implementation rerun after adopting `frscvNOMAD` `SIMPLE_LINE_SEARCH=yes` path default:

- `/tmp/crs_install_fr_simpleline_20260223.log`
- `/tmp/crs_nomad_post_frsimpleline_20260223_raw.csv`
- `/tmp/crs_nomad_post_frsimpleline_20260223_summary.csv`
- `/tmp/crs_nomad_post_frsimpleline_20260223_parity.rds`
- `/tmp/crs_nomad_post_frsimpleline_20260223.log`
- `/tmp/crs_nomad_frsimpleline_vs_pathdefaults_20260223_summary_compare.csv`
- `/tmp/crs_nomad_frsimpleline_vs_pre_20260223_summary_compare.csv`
- `/tmp/crs_nomad_exception_isolation_after_frdefaults_20260223_raw.csv`
- `/tmp/crs_nomad_exception_isolation_after_frdefaults_20260223_summary.csv`
- `/tmp/crs_nomad_exception_isolation_after_frdefaults_20260223.log`

## 2026-02-23 clean-break comparison vs NOMAD3 baseline

Reference pre-upgrade baseline:

- `/tmp/crs_nomad_pre_20260223_104449_raw.csv`
- `/tmp/crs_nomad_pre_20260223_104449_summary.csv`

Observed elapsed-time deltas (post vs pre, summary-level):

1. `frscvNOMAD`
   - fixed seed: mean `+665.85%`, median `+853.85%`
   - varying seed: mean `+954.84%`, median `+915.38%`
2. `krscvNOMAD`
   - fixed seed: mean `+578.64%`, median `+700.00%`
   - varying seed: mean `+981.82%`, median `+942.86%`
3. `npglpreg`
   - fixed seed: mean `-33.03%`, median `-34.36%`
   - varying seed: mean `+12.68%`, median `+33.23%`

Local-minima behavior across varying seeds (`npglpreg`, `frscvNOMAD`, `krscvNOMAD`):

1. `frscvNOMAD`: post objective better on `4/5` seeds
2. `krscvNOMAD`: post objective better on `3/5` seeds
3. `npglpreg`: post objective better on `5/5` seeds
   - best objective improved from `0.0421112191` (pre) to `0.036881` (post)
   - median objective improved from `0.0523705702` (pre) to `0.049461` (post)

Parity notes:

1. `frscvNOMAD` fixed-seed parameter signature changes: `0/5`
2. `krscvNOMAD` fixed-seed parameter signature changes: `0/5`
3. `npglpreg` fixed-seed parameter signature changes: `5/5` (drift remains)

Method note:

- `benchmarks/nomad/compare_prepost.R` can overcount fixed-seed repeats via cross-join.
- For this section, release-facing comparisons were taken from summary-level outputs and replicate-indexed parity/minima calculations.

## 2026-02-23 parameter-mapping reassessment

Question tested:

- Are the large kernel-path slowdowns mainly due to old-vs-new option mismatch?

Action:

1. Built an explicit NOMAD3.9.1-like compatibility profile in `snomadr()` using NOMAD4 names only.
2. Re-ran full benchmark suite and compared against both:
   - NOMAD3 baseline
   - prior clean-break NOMAD4 run without explicit profile

Result:

1. Relative to NOMAD3 baseline, kernel paths remained substantially slower:
   - `frscvNOMAD`: fixed mean roughly `+684%` to `+753%`, varying mean roughly `+985%` to `+1005%`
   - `krscvNOMAD`: fixed mean roughly `+630%` to `+719%`, varying mean roughly `+1230%` to `+1306%`
2. Relative to prior clean-break NOMAD4 run, this compatibility profile did not materially reduce runtimes:
   - mostly within low double-digit percent, often slightly slower
3. Local-minima behavior stayed broadly similar:
   - `npglpreg` still better on all varying-seed runs
   - `frscvNOMAD` mostly better
   - `krscvNOMAD` mixed

Conclusion:

- Parameter-name replacement (`MIN_POLL_SIZE` -> `MIN_FRAME_SIZE`) was not the main driver of observed runtime gaps.
- Remaining slowdown is likely dominated by NOMAD4 implementation/runtime characteristics on these small kernel-CV workloads rather than superficial option-name mismatch.

Refresh snapshot (`..._20260223_2...`):

1. Versus NOMAD3 baseline:
   - `frscvNOMAD`: fixed `+752.70%`, varying `+1005.49%` mean elapsed
   - `krscvNOMAD`: fixed `+719.21%`, varying `+1305.63%` mean elapsed
   - `npglpreg`: fixed `-25.38%`, varying `+19.52%` mean elapsed
2. Varying-seed local minima:
   - `frscvNOMAD`: post better `4/5`
   - `krscvNOMAD`: post better `2/5`
   - `npglpreg`: post better `5/5`

## 2026-02-23 path-specific default tuning outcome

Best robust profile found:

1. Keep global `snomadr()` compatibility profile as-is.
2. Add path defaults for `frscvNOMAD` and `krscvNOMAD` only:
   - `QUAD_MODEL_SEARCH = no`
   - `EVAL_QUEUE_SORT = DIR_LAST_SUCCESS`

Observed impact versus current NOMAD4 compatibility baseline:

1. `frscvNOMAD`: mean elapsed roughly `-56%` (fixed and varying), no objective drift observed.
2. `krscvNOMAD`: mean elapsed roughly `-55%` to `-57%`, objective drift negligible (about `2.4e-06` in varying-seed mean).
3. `npglpreg`: unchanged objective and essentially unchanged runtime (this path keeps existing MADS defaults).

Relative to NOMAD3 baseline after this tuning:

1. `frscvNOMAD`: slowdown reduced from roughly `+750% to +1005%` down to about `+292% to +358%`.
2. `krscvNOMAD`: slowdown reduced from roughly `+719% to +1306%` down to about `+240% to +498%`.
3. `npglpreg`: unchanged from prior NOMAD4 profile behavior.

Implemented-code snapshot (`...post_pathdefaults_20260223_1...`):

1. Versus prior NOMAD4 compatibility profile:
   - `frscvNOMAD`: about `-51%` to `-57%` mean elapsed
   - `krscvNOMAD`: about `-53%` to `-57%` mean elapsed
   - `npglpreg`: about `+2%` to `+4%` mean elapsed, objective unchanged in this run
2. Versus NOMAD3 baseline:
   - `frscvNOMAD`: about `+311%` (fixed), `+373%` (varying)
   - `krscvNOMAD`: about `+272%` (fixed), `+478%` (varying)
   - `npglpreg`: about `-22.8%` (fixed), `+23.3%` (varying)

## 2026-02-23 standalone solver-mode investigation

NOMAD4 provides many standalone optimization modes, but for `crs`/`npglpreg`:

1. Fast standalone modes (`CS_OPTIMIZATION`, `NM_OPTIMIZATION`, `QP_OPTIMIZATION`, `QUAD_MODEL_OPTIMIZATION`, `RANDOM_ALGO_OPTIMIZATION`) often produced materially worse objectives on at least one of the target paths (especially `npglpreg`) despite lower runtime.
2. OpenMP-gated modes (`PSD_MADS_OPTIMIZATION`, `COOP_MADS_OPTIMIZATION`) are unavailable in the current build.
3. Practical conclusion: MADS with path-specific search-option tuning remains the best current choice for balancing speed and solution quality stability.

## 2026-02-23 expanded solver/parameter space investigation

Goal:

- Test whether additional NOMAD4 solver families and tuned settings (`nmulti`, `max.bb.eval`, queue/search toggles) can beat current `crs` defaults on speed while preserving objective quality.

Artifacts:

1. Stage-1 broad scan:
   - `/tmp/crs_nomad_solver_space_stage1b_20260223_raw.csv`
   - `/tmp/crs_nomad_solver_space_stage1b_20260223_summary.csv`
   - `/tmp/crs_nomad_solver_space_stage1b_20260223_compare.csv`
   - `/tmp/crs_nomad_solver_space_stage1b_20260223_safe_rank.csv`
2. Stage-2 shortlist confirmation:
   - `/tmp/crs_nomad_solver_space_stage2_20260223_raw.csv`
   - `/tmp/crs_nomad_solver_space_stage2_20260223_compare.csv`
3. Stage-3 exploratory expansion (included exception-only profiles):
   - `/tmp/crs_nomad_solver_space_stage3_20260223_raw.csv`
   - `/tmp/crs_nomad_solver_space_stage3_20260223_compare.csv`
   - `/tmp/crs_nomad_solver_space_stage3_20260223_agg.csv`
4. Stage-3C clean expansion (full coverage, no exception profiles):
   - `/tmp/crs_nomad_solver_space_stage3c_20260223_raw.csv`
   - `/tmp/crs_nomad_solver_space_stage3c_20260223_compare.csv`
   - `/tmp/crs_nomad_solver_space_stage3c_20260223_agg.csv`
   - `/tmp/crs_nomad_solver_space_stage3c_20260223_strict_safe.csv`
   - `/tmp/crs_nomad_solver_space_stage3c_20260223.log`

Important method note:

1. Profiles that intentionally trigger NOMAD exceptions (for example unsupported `VNS_MADS_OPTIMIZATION` in this build, or DiscoMADS without revealing output) can contaminate multi-profile runs in a single R session by leaving conflicting algorithm flags.
2. Release-facing conclusions should therefore use clean-profile sweeps (Stage-3C) or isolate profiles per process when testing exception-prone modes.

Stage-3C strict-safe winners (`0` worsens, `0` errors, all `4` cells):

1. `frscvNOMAD`
   - `fr_current`: baseline
   - `fr_noquad_lexi`: mean elapsed `-2.49%`, worst cell `+4.17%`, objective unchanged
2. `krscvNOMAD`
   - `kr_current`: baseline
3. `npglpreg`
   - `np_current`: baseline
   - `np_quadbox1`: mixed speed (`+1.06%` mean; cell range about `-2.92%` to `+8.36%`) with net objective improvements (`5`) and no worsens

Stage-3C aggressive frontier (fast but objective risk):

1. `frscvNOMAD` with standalone CS/NM/QP/random:
   - elapsed roughly `-88%` to `-95%`
   - objective worsens present (`4/20` to `20/20` pairs depending on mode)
   - CS max mean objective drift observed about `+0.0019362`
2. `krscvNOMAD` with standalone CS/NM/QP/random:
   - elapsed roughly `-73%` to `-96%`
   - objective worsens present (`10/20` to `20/20`)
   - CS max mean objective drift observed about `+0.001487`
3. `npglpreg` aggressive profiles:
   - CS/NM/QP/random greatly faster (roughly `-39%` to `-98%`) but worsen almost always (`19/20` or `20/20`)
   - `np_noquad` gave useful speed (`-36%` mean) but still had mixed objective movement (`4` worsens, `7` improves)

Interpretation:

1. For production defaults, no broad solver-family switch improved speed while remaining strictly parity-safe across all tested cells.
2. Best robust default set remains:
   - `frscvNOMAD`: MADS with existing path defaults
   - `krscvNOMAD`: MADS with existing path defaults
   - `npglpreg`: current MADS compatibility profile (optionally evaluate `np_quadbox1` case-by-case)
3. Defaults currently in `crs` are still MADS-based; standalone algorithms are only activated by explicit options.

## 2026-02-23 extended solver-space follow-up and default update

Goal:

- Re-open the highest-value tuning tranche by testing additional MADS search-method controls while keeping strict parity/error gates.

Strict-safe outcomes from `/tmp/crs_nomad_solver_space_ext_20260223_strict_safe.csv`:

1. `frscvNOMAD`
   - `fr_simple_line`: mean elapsed `-5.06%`, worst cell `-0.017%`, `0` worsens, `0` errors
   - `fr_noquad_lexi`: mean elapsed `-4.27%`, worst cell `+4.71%`
   - `fr_vns_search`: mean elapsed `-2.49%`, worst cell `+15.77%`
   - `fr_current`: baseline
2. `krscvNOMAD`
   - `kr_current`: baseline remained best strict-safe profile
3. `npglpreg`
   - `np_current`: baseline remained best strict-safe profile
   - `np_quadbox1`: no worsens with `5` improves, but mean elapsed `+1.48%` (mixed speed effect)

Profiles rejected for production defaults:

1. Standalone solver modes (CS/NM/QP/random): very large speedups, but objective-worsen counts remained material.
2. `np_vns_search`: unstable in this build (`error_runs=19`).

Implemented code change from this follow-up:

1. `frscvNOMAD` path defaults now include:
   - `QUAD_MODEL_SEARCH=no`
   - `EVAL_QUEUE_SORT=DIR_LAST_SUCCESS`
   - `SIMPLE_LINE_SEARCH=yes`
   - `SPECULATIVE_SEARCH=no`
2. `krscvNOMAD` path defaults remain:
   - `QUAD_MODEL_SEARCH=no`
   - `EVAL_QUEUE_SORT=DIR_LAST_SUCCESS`

Post-update checkpoint:

1. Versus prior path-default build (`/tmp/crs_nomad_frsimpleline_vs_pathdefaults_20260223_summary_compare.csv`):
   - `frscvNOMAD`: fixed `-4.77%`, varying `+0.24%`, objective/parameter diffs `0`
2. Versus NOMAD3 baseline (`/tmp/crs_nomad_frsimpleline_vs_pre_20260223_summary_compare.csv`):
   - `frscvNOMAD`: fixed `+291.60%`, varying `+373.38%`
   - `krscvNOMAD`: fixed `+263.32%`, varying `+463.08%`
3. Exception-isolation recheck after default change remained clean:
   - `cs_status_in = ok`, `nm_status_in = ok` across tested scenario/seed pairs.

## 2026-02-23 MADS-only deep sweep (cross-`n` confirmation)

Goal:

- Expand tuning beyond search toggles and test direction/restart/budget controls under strict parity/error gates, then confirm winners at a larger sample size.

Method:

1. Run `run_nomad_mads_deep_space.R` at `n=180` and `n=300`.
2. Use strict-safe filter (`cells_present=4`, `worsen_count=0`, `error_runs=0`).
3. Keep only profiles that remain strict-safe across both runs.

Cross-`n` strict-safe outcomes:

1. `frscvNOMAD`
   - `fr_dirnp1neg` remained best:
     - `n=180`: mean elapsed `-82.04%`
     - `n=300`: mean elapsed `-81.20%`
     - objective: no worsens, no improves in all four cells
   - `fr_dir2n` was close second and also strict-safe.
2. `krscvNOMAD`
   - `kr_nmtrial40` remained best strict-safe:
     - `n=180`: mean elapsed `-4.32%`
     - `n=300`: mean elapsed `-2.52%`
     - objective: no worsens in all four cells
3. `npglpreg`
   - no robust speed gain over baseline survived cross-`n` confirmation.
   - `np_nmtrial40` remained strict-safe but effectively neutral (`~0.04%` mean at `n=300`).

Implemented default update from this sweep:

1. `frscvNOMAD` path defaults now additionally set:
   - `DIRECTION_TYPE=ORTHO N+1 NEG`
2. `krscvNOMAD` defaults left unchanged after canonical-bench validation:
   - `NM_SEARCH_MAX_TRIAL_PTS_NFACTOR=40` remains a validated optional profile (faster in deep mixed/hard sweeps, mixed in canonical bench).
3. `npglpreg` defaults unchanged.

Implementation note:

1. Profiles using `SIMPLE_LINE_SEARCH=yes` with speculative search still enabled are invalid in NOMAD4 and triggered expected exceptions:
   - `SimpleLineSearchMethod: cannot work with speculative search`
2. These invalid profiles are retained as rejected candidates in deep-sweep artifacts and excluded by strict-safe filters.

## 2026-02-23 exception-path contamination investigation and fix

Question tested:

- Do exception-only NOMAD modes contaminate later optimizer calls in the same R session?

Pre-fix evidence (before C-interface cleanup patch):

1. Study log:
   - `/tmp/crs_nomad_exception_isolation_study_20260223.log`
2. Observed behavior:
   - after a `DISCO_MADS_OPTIMIZATION` failure, subsequent in-session `npglpreg` runs (`cs_opt`, `nm_opt`, and even baseline) failed,
   - matching isolated subprocess calls remained `ok`.
3. Summary signature from that run:
   - `cs_status_in = error`, `cs_status_isolated = ok` for all tested scenario/seed pairs
   - `nm_status_in = error`, `nm_status_isolated = ok` for all tested scenario/seed pairs

Root cause:

1. `solveNomadProblem()` in
   - `/Users/jracine/Development/crs/src/nomad4_src/interfaces/CInterface/NomadStdCInterface.cpp`
2. success path reset global NOMAD components (`MainStep::resetComponentsBetweenOptimization`) but exception paths did not, leaving algorithm-state contamination for later calls.

Fix implemented:

1. Added shared cleanup routine in `solveNomadProblem()` and invoked it on all return paths (including `checkAndComply` failures and run-time exceptions).
2. Cleanup made exception-safe to avoid failures when cache is not yet instantiated.

Post-fix evidence:

1. Study artifacts:
   - `/tmp/crs_nomad_exception_isolation_study_20260223_afterfix_run.log`
   - `/tmp/crs_nomad_exception_isolation_study_20260223_afterfix_raw.csv`
   - `/tmp/crs_nomad_exception_isolation_study_20260223_afterfix_summary.csv`
2. Outcome:
   - `cs_status_in = ok` and `nm_status_in = ok` for all tested scenario/seed pairs after a `disco_opt` failure,
   - in-session baseline calls recovered and remained runnable.
3. Regression smoke after patch:
   - `/tmp/crs_nomad_smoke_after_exception_reset_20260223_summary.csv`
   - `/tmp/crs_nomad_smoke_after_exception_reset_20260223_parity.rds`

## 2026-02-23 OpenMP/parallel-evaluation investigation (MacPorts `libomp`)

Question tested:

- Can NOMAD4 OpenMP parallel evaluation (`NB_THREADS_PARALLEL_EVAL`) be integrated in `crs` and materially speed optimization?

Build/integration artifacts:

1. OpenMP-enabled install (experimental local build):
   - `/tmp/crs_install_openmp_exp_20260223_retry3.log`
2. Clean non-OpenMP rebuild (true baseline after object purge):
   - `/tmp/crs_install_noopenmp_clean_20260223.log`

Runtime safety probe artifacts:

1. Per-case thread matrix logs:
   - `/tmp/crs_openmp_probe_fr_t*.log`
   - `/tmp/crs_openmp_probe_kr_t*.log`
   - `/tmp/crs_openmp_probe_np_t*.log`
   - `/tmp/crs_openmp_probe_basic_t*.log`
2. Matrix summaries:
   - `/tmp/crs_openmp_probe_summary.psv`
   - `/tmp/crs_openmp_probe_analysis.csv`

OpenMP build-overhead artifacts (`threads=1` only):

1. OpenMP build benchmark:
   - `/tmp/crs_nomad_openmp_t1_20260223_raw.csv`
   - `/tmp/crs_nomad_openmp_t1_20260223_summary.csv`
   - `/tmp/crs_nomad_openmp_t1_20260223.log`
2. Clean non-OpenMP benchmark:
   - `/tmp/crs_nomad_noopenmp_clean_20260223_raw.csv`
   - `/tmp/crs_nomad_noopenmp_clean_20260223_summary.csv`
   - `/tmp/crs_nomad_noopenmp_clean_20260223.log`

Findings:

1. Integration is technically feasible for local experimentation (OpenMP compile/link works), but current `crs` callback architecture is not safe for `NB_THREADS_PARALLEL_EVAL > 1`.
2. Safety matrix outcome:
   - `threads=1`: stable for `frscvNOMAD`, `krscvNOMAD`, `npglpreg`, `snomadr` basic case
   - `threads>=2`: mostly fatal (`SIGSEGV`, `SIGBUS`, abort) or corrupted behavior
   - one apparent `fr` "`ok`" run at high thread count returned `objective=0` while log contained internal R errors; treated as invalid/corrupted result
3. Performance impact of OpenMP build alone (still `threads=1`) was negative:
   - about `+5.25%` mean elapsed overall versus clean non-OpenMP build on canonical harness
   - objective parity remained exact in this comparison

Decision:

1. Keep production `crs` on non-OpenMP NOMAD build for now.
2. Do not recommend or document `NB_THREADS_PARALLEL_EVAL > 1` for current `crs`/`npglpreg` paths.
3. Revisit only if evaluation callbacks are moved to a thread-safe native path that avoids concurrent R API re-entry.

## 2026-02-23 `krscvNOMAD` time-budget frontier (strict vs aggressive)

Goal:

- Determine whether a much faster `krscvNOMAD` profile can improve practical global-search outcomes, and whether any such profile is robust enough for defaults.

Artifacts:

1. Frontier sweep at `n=300`:
   - `/tmp/crs_nomad_kr_timebudget_20260223_n300_raw.csv`
   - `/tmp/crs_nomad_kr_timebudget_20260223_n300_compare.csv`
   - `/tmp/crs_nomad_kr_timebudget_20260223_n300_agg.csv`
   - `/tmp/crs_nomad_kr_timebudget_20260223_n300_strict_safe.csv`
   - `/tmp/crs_nomad_kr_timebudget_20260223_n300.log`
2. Frontier sweep at `n=500`:
   - `/tmp/crs_nomad_kr_timebudget_20260223_n500_raw.csv`
   - `/tmp/crs_nomad_kr_timebudget_20260223_n500_compare.csv`
   - `/tmp/crs_nomad_kr_timebudget_20260223_n500_agg.csv`
   - `/tmp/crs_nomad_kr_timebudget_20260223_n500_strict_safe.csv`
   - `/tmp/crs_nomad_kr_timebudget_20260223_n500.log`
3. Targeted stress runs:
   - `/tmp/crs_kr_targeted_n500_raw.csv`
   - `/tmp/crs_kr_targeted_n500_compare.csv`
   - `/tmp/crs_kr_dir2n_eval140_vs_current_n300_50seeds.csv`

Key profiles tested:

1. Baseline:
   - `kr_current` (`nmulti=2`, `max.bb.eval=100`, current defaults)
2. Strict-safe candidate:
   - `kr_nmtrial40` (`opts = list(NM_SEARCH_MAX_TRIAL_PTS_NFACTOR=40)`)
3. Aggressive candidate:
   - `kr_dir2n_eval140` (`opts = list(DIRECTION_TYPE="ORTHO 2N")`, `max.bb.eval=140`)

Findings:

1. `kr_nmtrial40` stayed strict-safe across tested cells with small speed gain (roughly `-0.6%` to `-1.3%` mean).
2. `kr_dir2n_eval140` delivered large speed gains:
   - around `-58.7%` mean (`n=300` frontier)
   - around `-57.0%` mean (`n=500` frontier)
3. Objective behavior for `kr_dir2n_eval140`:
   - often improved, especially on hard scenarios,
   - but not strictly parity-safe in larger-seed stress runs:
     - `n=300`, 50-seed hard scenario: `26` improves, `10` worsens, `14` ties (`tol=1e-6`)
     - worst observed positive drift about `+2.4e-05`, best improvement about `-1.432e-03`
4. Interpretation:
   - this profile is a strong aggressive search mode (fast and often better minima),
   - but not suitable as a strict default under current parity gates.
5. Additional robustness note:
   - `kr_simple_specoff` showed sporadic runtime exceptions in targeted stress runs (`SimpleLineSearchMethod: evaluated point not found in cache`), so it is not a reliable candidate despite occasional speed gains.

Decision:

1. Keep `krscvNOMAD` package defaults conservative.
2. Keep `kr_nmtrial40` as optional strict-safe speed tweak.
3. Keep `kr_dir2n_eval140` documented as an aggressive, user-opt-in exploration profile for faster global-search attempts.

## Current gate state

1. Build/install: passing.
2. Runtime compatibility: passing for core paths (`frscvNOMAD`, `krscvNOMAD`, `npglpreg`, `snomadr`).
3. Strict parity/performance gates: not yet fully passing.

Outstanding gaps:

1. Fixed-seed performance regressions remain in `frscvNOMAD` and `krscvNOMAD` vs NOMAD3 baseline.
2. `npglpreg` still shows parameter drift under strict comparisons.
3. Direct `snomadr` fixed-seed behavior is sensitive to search-profile changes.

## Practical forward plan

1. Keep bridge contract stable and continue tuning at option-profile level.
2. Prefer path-aware tuning over global search-disable flags.
3. Re-run strict pre/post gates after each tuning checkpoint and append artifacts here.
4. Treat this file as the only canonical place for migration status and gate outcomes.
