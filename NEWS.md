# crs 0.15-44

* Exposed NOMAD native cache counters from `snomadr()` (`cache.hits`,
  `cache.size`, `callback.evaluations`, and `total.evaluations`) so downstream
  callers can distinguish true R callback evaluations from points satisfied by
  NOMAD's cache without changing optimizer behavior.
