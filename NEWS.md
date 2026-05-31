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
