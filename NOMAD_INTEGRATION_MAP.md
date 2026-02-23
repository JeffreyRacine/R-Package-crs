# NOMAD Integration Map (`crs`)

Snapshot date: 2026-02-23

This file is intentionally static/technical.

Dynamic run outcomes and gate status live in:

- `/Users/jracine/Development/crs/NOMAD_STATUS_AND_GATES.md`

## Current embedded NOMAD version (wired in `crs`)

`4.5.0`

Evidence in source:

- `/Users/jracine/Development/crs/src/nomad4_src/src/nomad_version.hpp` defines:
  - `NOMAD` namespace alias as `NOMAD_4_5`
  - `NOMAD_VERSION_NUMBER "4.5.0"`
- `/Users/jracine/Development/crs/src/Makevars` compiles from `./nomad4_src`
- `/Users/jracine/Development/crs/src/snomadr.cpp` includes:
  - `NomadStdCInterface.h`
  - `nomad_version.hpp`

## Clean-break note

Legacy NOMAD3 source tree `src/nomad_src` has been removed from `crs`.
The package now keeps only the NOMAD4 codebase and C interface bridge.

## Integration topology (current)

1. Build-time vendoring and linkage

- Embedded NOMAD4 source is compiled directly from:
  - `/Users/jracine/Development/crs/src/nomad4_src/src`
  - `/Users/jracine/Development/crs/src/nomad4_src/ext/sgtelib/src`
  - `/Users/jracine/Development/crs/src/nomad4_src/interfaces/CInterface/NomadStdCInterface.cpp`
- Build graph is defined in:
  - `/Users/jracine/Development/crs/src/Makevars`

2. Native routine registration into R

- `/Users/jracine/Development/crs/NAMESPACE` loads native routines via `useDynLib(crs, .registration = TRUE)`
- `/Users/jracine/Development/crs/src/crs_init.c` registers:
  - `snomadRInfo`
  - `snomadRSolve`
  - `smultinomadRSolve`

3. R-level solver wrapper API

- Main wrapper:
  - `/Users/jracine/Development/crs/R/snomadr.R`
- Calls into native bridge:
  - `.Call(snomadRInfo, ...)`
  - `.Call(snomadRSolve, ...)`
  - `.Call(smultinomadRSolve, ...)`

4. C++ bridge and callback path

- Bridge implementation:
  - `/Users/jracine/Development/crs/src/snomadr.cpp`
- Uses NOMAD4 C interface:
  - `createNomadProblem(...)`
  - `addNomad*Param(...)`
  - `solveNomadProblem(...)`
  - `createNomadResult(...)` / `load*SolutionsNomadResult(...)`
- Objective callback re-enters R via `R_tryEval(...)`.

5. High-level `crs` entry points that route into NOMAD

- `frscvNOMAD(...)`
- `krscvNOMAD(...)`
- `npglpreg(...)` (via NOMAD-driven CV path)
- additional direct `snomadr(...)` callers in package code.

## Current bridge behavior for options

- `snomadr(opts=...)` passes option names/values to NOMAD4.
- `crs` expects NOMAD4 option names directly (for example `MIN_FRAME_SIZE`).
- array option normalization for NOMAD4 format remains in place.
- relative (`r...`) scalar parsing for array options remains in place.
- integer/granular mesh/frame floor handling remains in place.
- `MAX_EVAL` auto-set from `MAX_BB_EVAL` when absent.
- safe `EPSILON` sanitization for NOMAD4 minimum precision constraints.

Canonical option catalog for embedded NOMAD4:

- `/Users/jracine/Development/crs/inst/nomad/NOMAD_4_5_0_OPTIONS_REFERENCE.md`

## Working call flow

1. R caller builds NOMAD problem state (`n`, `bbin`, `bbout`, bounds, `x0`, options, callback env).
2. `snomadr(...)` dispatches to `.Call(...)`.
3. C++ bridge constructs and parameterizes a NOMAD4 C-interface problem.
4. NOMAD4 evaluates trial points through the callback into R objective code.
5. Bridge maps solver outputs/status to the legacy R-facing return contract:
   - `status`, `message`, `bbe`, `iterations`, `objective`, `solution`.

## Practical implications

- `crs` is now wired to embedded NOMAD4 (`4.5.0`), not NOMAD3 (`3.9.1`).
- `crs` keeps a clean-break NOMAD4-only source layout.
