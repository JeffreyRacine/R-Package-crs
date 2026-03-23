# NOMAD Markdown Consolidation Closeout (2026-03-22)

Purpose:

- reduce overlapping root-level NOMAD notes in `crs`,
- keep one live canonical NOMAD note in the repo root,
- preserve superseded companion notes in a dated archive bundle.

Live canonical note after consolidation:

- `/Users/jracine/Development/crs/NOMAD_STATUS_AND_GATES.md`

Files archived here:

1. `NOMAD_FORWARD_GUIDE.md`
   - superseded as a separate root document because its active operational
     guidance is now folded into the canonical status file.
2. `NOMAD_INTEGRATION_MAP.md`
   - superseded as a separate root document because the still-relevant
     integration snapshot and OpenMP implications are now summarized in the
     canonical status file.

Files intentionally kept live in the repo root:

1. `/Users/jracine/Development/crs/NOMAD_STATUS_AND_GATES.md`
   - active canonical NOMAD architecture/status/gates/OpenMP note
2. `/Users/jracine/Development/crs/NOMAD_391_TO_450_OPTION_MAP.md`
   - durable option/provenance mapping reference
3. `/Users/jracine/Development/crs/CRS_MODERNIZATION_PLAN.md`
   - broader package modernization tracker

Follow-up captured in the canonical note:

1. current experimental OpenMP path can remain in-tree without requiring a fork
2. ordinary `snomadr(eval.f=...)` remains unsafe for threaded NOMAD evaluation
   because callbacks re-enter R
3. the same limitation currently carries over to `np` and `npRmpi` routes that
   call `crs::snomadr()` with wrapped R closures
