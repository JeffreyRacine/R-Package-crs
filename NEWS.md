# crs 0.15-46

* Modernized `crsiv()` and `crsivderiv()` with training-grid formula/data
  interfaces (`y ~ z | w` and `y ~ z | w | x`), one authoritative row map for
  subset/NA/weights/starting values, collision-safe internal role names, and
  namespaced IV metadata. Added selected-state `fitted()`, `residuals()`, and
  `predict()` semantics, structured summaries, corrected original-response
  plot overlays, and self-contained formula-based examples. Native evaluation
  interfaces and established CRS post-fit projection remain available;
  formula-time evaluation and formula-object `newdata` are deliberately
  deferred and fail clearly.

* Corrected fitted-value centering throughout `crsivderiv()` so both empirical
  terms of the Equation (14) adjoint use the same fitted conditional-residual
  vector. The initial adjoint had centered its second term on the raw
  residual, while the constructed-residual route also averaged the wrong
  constructed quantity in both initial and recursive states. This can
  materially change derivative trajectories, stopping states, and fitted
  curves, and removes catastrophic finite-sample losses observed under the
  malformed constructed-residual adjoint.

* Made `crsivderiv()` Landweber-Fridman states coherent with the estimator
  definition: iteration `N` now stores the matched derivative, integrated
  curve, and stopping criterion after exactly `N` updates. Explicit
  `starting.values` are the state used immediately before update 1, without a
  hidden preliminary update. Returned `phi` and `phi.prime`, the reported
  iteration, and the summary stopping value now all refer to the selected
  state. This corrects the previous one-column mismatch between the curve and
  derivative histories; results can change when that mismatch affected state
  selection or extraction.

* Corrected the private Gaussian integral operator used by `crsivderiv()` so
  the empirical adjoint applies the ordinary kernel CDF required by Equation
  (14) of Florens, Centorrino, and Racine. The helper had multiplied the CDF
  average by the continuous bandwidth, producing the same malformed scaling
  implicated in `np` issue #57. Other spline fitting, bandwidth, stopping,
  optimizer, and public helper behavior is unchanged.

* Finalized the observed native NOMAD interface so explicit user interruption
  has a dedicated status after orderly native cleanup, ordinary observer
  errors remain fail-open, and a failed callback cannot expose unwritten
  black-box output memory. Existing native ABI values, layouts, callable names,
  and successful-run numerical/accounting behavior are unchanged.

* Hardened native NOMAD interruption by limiting the operating-system signal
  handler to signal-safe state capture, consuming interruption from ordinary
  solver control flow, and restoring the prior process handler before return.
  Interrupt status now takes precedence when a solve is interrupted after an
  earlier callback failure, and R-callback failures retain their first trapped
  diagnostic in the bounded result message.

# crs 0.15-45

* Improved mixed-data cross-validation efficiency by using a guarded
  weighted least-squares Gram/Cholesky solve with QR/SVD fallback for
  kernel and factor spline CV routes. The change covers `cv.ls`, `cv.gcv`,
  and `cv.aic` for additive, tensor, and GLP bases consistently across
  weighted and unweighted mean-regression CV, while preserving the existing
  rank, fallback, and objective contracts.

* Restored the GLP model-matrix column oracle used to select generalized
  local-polynomial interactions while retaining compiled column-product
  construction, preserving GLP stability and shape compatibility.

* Improved categorical-kernel mean fit and evaluation efficiency by reusing
  the guarded weighted least-squares Gram/Cholesky primitive for the
  non-quantile prediction ingredients when `model.return=FALSE`, covering
  additive, tensor, and GLP bases consistently across weighted and
  unweighted routes while preserving fitted values, intervals, standard
  errors, hatvalues, rank, and residual degrees of freedom.

* Improved categorical-kernel mean-gradient efficiency by reusing the same
  guarded weighted least-squares Gram/Cholesky primitive for non-quantile
  derivative solves in training and evaluation routes. The change covers
  additive, tensor, and GLP bases consistently across weighted and unweighted
  routes while preserving derivative estimates and normal-approximation
  intervals to numerical roundoff.

* Changed the public `crs.formula()` NOMAD evaluation-budget default
  `max.bb.eval` to `NULL`, allowing `crs()` to choose route-specific
  defaults after route selection. Continuous-only `frscvNOMAD` searches
  now default to `max.bb.eval=10000`, while kernel/categorical
  `krscvNOMAD` searches default to `max.bb.eval=1000`. These defaults
  were set on the basis of simulation evidence and real-world
  applications; explicit user-supplied values continue to override the
  route defaults.

* Repaired public-wrapper consistency so `crs()` forwards NOMAD search
  controls to `frscvNOMAD` consistently with `krscvNOMAD`, including
  `random.seed`, `max.bb.eval`, integer mesh/frame geometry controls,
  and quantile level `tau`.

* Made NOMAD multistart generation invariant to progress/display
  settings in `frscvNOMAD` and `krscvNOMAD`, so enabling or suppressing
  optimizer progress no longer changes the starting-value geometry or
  fitted result.

* Implemented independent native NOMAD restart sweeps for CRS
  cross-validation searches, preserving the public multistart contract while
  keeping restart diagnostics explicit.

* Added public `max.eval` control for NOMAD point-lookup budgets and changed
  continuous-only `frscvNOMAD` defaults so `MAX_EVAL` and `MAX_BB_EVAL` can be
  controlled independently. The default point-lookup cap for continuous-only
  CRS NOMAD searches is now 1000, reducing duplicate NOMAD cache lookups while
  preserving the established black-box evaluation budget.

* Clarified NOMAD cache reporting in summaries. User-facing output now
  distinguishes true objective/function evaluations from repeated NOMAD point
  lookups avoided by the native cache.

* Added elapsed-time recording and reporting for CRS cross-validation and
  summary paths using wall-clock elapsed time rather than user CPU time.

* Modernized CRS plot methods toward the current `np` plot interface.
  `plot.crs()` now displays fitted mean/quantile functions by default,
  accepts NP-style controls such as `errors`, `band`, `B`, `output`,
  `data_overlay`, `data_rug`, `perspective`, `renderer`, and `view`, and
  supports transparent viridis base surfaces, base rotation, rgl surface
  extras, data overlays, rugs, and fitted-surface asymptotic/inid-bootstrap
  intervals. Legacy CRS plot switches such as `ci`, `mean`, `plot.rug`,
  and `plot.errors.*` now fail fast with NP-interface guidance.

* Hardened CRS plot argument validation and public plot contracts. Unsupported
  controls now fail explicitly rather than being silently ignored, categorical
  plot rendering and legends are aligned more closely with `np`, and CRS
  gradient plots now support refit and wild-bootstrap intervals where defined.

* Repaired explicit derivative handling in `predict.crs()` so caller-supplied
  derivative orders are honored rather than bypassed by the fitted-object fast
  path.

* Expanded the `crs_nomad_api` help page with package-author guidance,
  direct C-callback and R-callback bridge skeletons, and explanations of the
  most important native result fields.

# crs 0.15-44

* Added the final package-author native NOMAD C API
  (`crs_nomad_solve`) with C and R callback modes, explicit option arrays,
  generated starts, multistart support, budget handling, mesh-size array
  options, interrupt handling, and deterministic seed behavior for downstream
  packages.
* Exposed NOMAD native cache counters from `snomadr()` (`cache.hits`,
  `cache.size`, `callback.evaluations`, and `total.evaluations`) so downstream
  callers can distinguish true R callback evaluations from points satisfied by
  NOMAD's cache without changing optimizer behavior.
* Hardened native NOMAD API contracts by rejecting cache-off native solves,
  mirroring infinite objective handling, documenting main-thread affinity for
  R callback routes, and tightening named-list validation in native helpers.
