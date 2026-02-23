# NOMAD Integration Map (`crs`)

Snapshot date: 2026-02-23

## Current embedded NOMAD version (wired in `crs`)

`3.9.1`

Evidence in source:

- `/Users/jracine/Development/crs/src/nomad_src/nomad_version.hpp:45` sets `MAJOR_VER 3`
- `/Users/jracine/Development/crs/src/nomad_src/nomad_version.hpp:46` sets `MINOR_VER 9`
- `/Users/jracine/Development/crs/src/nomad_src/nomad_version.hpp:47` sets `REV_VER 1`
- `/Users/jracine/Development/crs/src/nomad_src/nomad_version.hpp:57` defines `NOMAD_VERSION_NUMBER` from those macros
- `/Users/jracine/Development/crs/src/nomad_src/defines.hpp:155` uses `NOMAD_VERSION_NUMBER` as `BASE_VERSION`
- `/Users/jracine/Development/crs/src/nomad_src/nomad.cpp:335` prints `NOMAD::VERSION`

Binary confirmation:

- `strings /Users/jracine/Development/crs/src/crs.so` contains `3.9.1` and `NOMAD version 3`.

## Integration topology

1. Build-time vendoring and linkage

- Embedded sources are compiled directly from:
  - `/Users/jracine/Development/crs/src/nomad_src`
  - `/Users/jracine/Development/crs/src/sgtelib_src`
- `/Users/jracine/Development/crs/src/Makevars:1` defines `NOMAD_SRC`
- `/Users/jracine/Development/crs/src/Makevars:2` defines `SGTELIB_SRC`
- `/Users/jracine/Development/crs/src/Makevars:3` onward adds NOMAD/SGTELIB objects to `OBJECTS`
- `/Users/jracine/Development/crs/src/Makevars:24` adds include paths for both trees

2. Native routine registration into R

- `/Users/jracine/Development/crs/NAMESPACE:1` loads shared library via `useDynLib(crs, .registration = TRUE)`
- `/Users/jracine/Development/crs/src/crs_init.c:17` declares `smultinomadRSolve`
- `/Users/jracine/Development/crs/src/crs_init.c:18` declares `snomadRInfo`
- `/Users/jracine/Development/crs/src/crs_init.c:19` declares `snomadRSolve`
- `/Users/jracine/Development/crs/src/crs_init.c:32` registers `smultinomadRSolve`
- `/Users/jracine/Development/crs/src/crs_init.c:33` registers `snomadRInfo`
- `/Users/jracine/Development/crs/src/crs_init.c:34` registers `snomadRSolve`

3. R-level solver wrapper API

- Primary R wrapper:
  - `/Users/jracine/Development/crs/R/snomadr.R:42`
- Info/version/help path:
  - `/Users/jracine/Development/crs/R/snomadr.R:85` calls `.Call(snomadRInfo, ...)`
- Solve path:
  - `/Users/jracine/Development/crs/R/snomadr.R:192` calls `.Call(snomadRSolve, ...)` for single run
  - `/Users/jracine/Development/crs/R/snomadr.R:194` calls `.Call(smultinomadRSolve, ...)` for multi-start
- Option marshaling helper:
  - `/Users/jracine/Development/crs/R/get.option.types.R:32`

4. C++ bridge and objective callback

- Bridge implementation:
  - `/Users/jracine/Development/crs/src/snomadr.cpp`
- R objective is evaluated inside C++ evaluator:
  - `/Users/jracine/Development/crs/src/snomadr.cpp:180` (`eval_f`)
  - `/Users/jracine/Development/crs/src/snomadr.cpp:218` (`RMy_Evaluator::eval_x`)
- R option list is converted to NOMAD parameters:
  - `/Users/jracine/Development/crs/src/snomadr.cpp:77` (`setApplicationOptions`)
  - `/Users/jracine/Development/crs/src/snomadr.cpp:140` uses `p.read(fp)`
- Runtime `nomad.opt` override is supported:
  - `/Users/jracine/Development/crs/src/snomadr.cpp:491`
  - `/Users/jracine/Development/crs/src/snomadr.cpp:764`

5. High-level `crs` entry points that route into NOMAD

- Core spline CV dispatch:
  - `/Users/jracine/Development/crs/R/crs.R:560` calls `frscvNOMAD(...)` when `cv == "nomad"`
  - `/Users/jracine/Development/crs/R/crs.R:620` calls `krscvNOMAD(...)` when `cv == "nomad"`
- NOMAD CV wrappers call `snomadr(...)`:
  - `/Users/jracine/Development/crs/R/frscvNOMAD.R:307`
  - `/Users/jracine/Development/crs/R/krscvNOMAD.R:370`
- Additional direct NOMAD usage:
  - `/Users/jracine/Development/crs/R/clsd.R:836`
  - `/Users/jracine/Development/crs/R/np.regression.glp.R:2539`
  - `/Users/jracine/Development/crs/R/np.regression.glp.R:2568`

## Working call flow

1. User calls `crs(...)` (or `glpcvNOMAD(...)`, `clsd(..., NOMAD=TRUE)`, etc.)
2. R code builds NOMAD search problem (`x0`, bounds, variable types, options, objective wrapper)
3. `snomadr(...)` sends packed state to `.Call(snomadRSolve)` or `.Call(smultinomadRSolve)`
4. C++ bridge creates `NOMAD::Parameters`, applies options, and runs `NOMAD::Mads`
5. For each trial point, C++ callback re-enters R objective (`eval.f`)
6. C++ returns status/objective/solution to R as list
7. R caller post-processes solution into model/CV outputs

## Practical implications of current state

- `crs` is currently wired to the legacy NOMAD v3 codebase (3.9.1), not NOMAD v4.x.
- Integration is source-vendored, not dynamically linked to a system NOMAD.
- Updating NOMAD in `crs` requires both:
  - native bridge compatibility work (`snomadr.cpp` / API differences), and
  - build graph updates in `/Users/jracine/Development/crs/src/Makevars`.
