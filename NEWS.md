# crs 0.15-45

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
