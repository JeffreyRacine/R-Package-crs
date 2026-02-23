# NOMAD 3.9.1 -> 4.5.0 Option Map (`crs`)

Last updated: 2026-02-23

Purpose:

- Provide a practical mapping from legacy NOMAD3 parameter names to NOMAD4 names.
- Record which NOMAD3 defaults are explicitly emulated in `crs`.

Primary references used:

- NOMAD 3.9.1 User Guide (PDF): https://www.gerad.ca/software/nomad/Downloads/user_guide.pdf
- Legacy NOMAD3 source defaults (historical `src/nomad_src/Parameters.cpp` in this repo history)
- Embedded NOMAD4 definitions: `src/nomad4_src/src/Attribute/*AttributesDefinition*.txt`

## Key mappings

| NOMAD 3.9.1 | NOMAD 4.5.0 | Notes |
| --- | --- | --- |
| `INITIAL_POLL_SIZE` | `INITIAL_FRAME_SIZE` | Direct conceptual replacement in NOMAD4. |
| `MIN_POLL_SIZE` | `MIN_FRAME_SIZE` | Direct conceptual replacement in NOMAD4. |
| `MODEL_SEARCH` | `QUAD_MODEL_SEARCH` + `SGTELIB_MODEL_SEARCH` | `MODEL_SEARCH QUADRATIC` corresponds to `QUAD_MODEL_SEARCH yes`, `SGTELIB_MODEL_SEARCH no`. |
| `MODEL_EVAL_SORT` | `EVAL_QUEUE_SORT` | `MODEL_EVAL_SORT quadratic` corresponds to `EVAL_QUEUE_SORT QUADRATIC_MODEL`. |
| `OPPORTUNISTIC_EVAL` | `EVAL_OPPORTUNISTIC` | Same high-level meaning. |
| `MODEL_RADIUS_FACTOR` | `QUAD_MODEL_SEARCH_BOX_FACTOR` and `QUAD_MODEL_BOX_FACTOR` | NOMAD4 split behavior into search/model box factors. |
| `MAX_EVAL_INTENSIFICATION` | `QUAD_MODEL_MAX_EVAL` + `SGTELIB_MODEL_MAX_EVAL` | NOMAD4 replacement per deprecated definition file. |

## Items with no strict one-to-one behavior

1. `MODEL_SEARCH_MAX_TRIAL_PTS`: marked as not implemented in NOMAD4 deprecated map.
2. `MODEL_SEARCH_OPTIMISTIC` / `MODEL_SEARCH_PROJ_TO_MESH`: marked as not implemented in NOMAD4 deprecated map.
3. `MESH_TYPE`: marked deprecated in NOMAD4 (internal mesh architecture changed).

## `crs` compatibility profile (NOMAD4 names only)

`snomadr()` now injects the following defaults if user `opts` does not set them:

1. `QUAD_MODEL_SEARCH = yes`
2. `SGTELIB_MODEL_SEARCH = no`
3. `NM_SEARCH = yes`
4. `SPECULATIVE_SEARCH = yes`
5. `EVAL_OPPORTUNISTIC = yes`
6. `EVAL_QUEUE_SORT = QUADRATIC_MODEL`
7. `DIRECTION_TYPE = ORTHO N+1 QUAD`
8. `QUAD_MODEL_SEARCH_BOX_FACTOR = 2.0`
9. `QUAD_MODEL_BOX_FACTOR = 2.0`

Design rule:

- User-supplied `opts` always override these defaults.
- `crs` uses NOMAD4 names directly; no legacy-name alias translation is reintroduced.

## Optimizer continuity: NOMAD3 vs NOMAD4

1. Both versions are from the same NOMAD project and remain MADS-family derivative-free black-box optimizers.
2. NOMAD4 is not a drop-in behavior clone of NOMAD3.9.1:
   - internal implementation and runtime architecture changed,
   - several NOMAD3 options are deprecated or not implemented in NOMAD4,
   - some model/search controls were split into multiple NOMAD4 options.
3. Practical implication for `crs`: use NOMAD4 option names and benchmark gates for behavior validation, not assumption of strict NOMAD3 runtime equivalence.
