# crs Rd Argument And Documentation Organization Campaign Plan

Date: 2026-04-28

## Goal

Make the `crs` package documentation easier for users to navigate by organizing
function arguments by user mental model, clarifying pass-through arguments, and
checking that `.Rd` files faithfully describe the current implementation, while
preserving the public API and numerical behavior.

This is a documentation-only campaign. Actual R function formal reordering is
not needed for the goal and is out of scope. The `.Rd` files should be made as
orderly and useful as possible while preserving the package API exactly.

## Scope

Primary documentation targets:

- `man/crs.Rd`
- `man/clsd.Rd`
- `man/snomadr.Rd`
- `man/frscv.Rd`
- `man/frscvNOMAD.Rd`
- `man/krscv.Rd`
- `man/krscvNOMAD.Rd`
- `man/crsiv.Rd`
- `man/crsivderiv.Rd`
- `man/crssigtest.Rd`

Secondary documentation targets:

- `man/gsl-bs.Rd`
- `man/glp.model.matrix.Rd`
- `man/tensor.prod.model.matrix.Rd`
- `man/uniquecombs.Rd`
- package and data pages, only if light cleanup is clearly helpful

Out of scope for this tranche:

- changing estimator behavior,
- changing default argument values,
- changing C/C++ or NOMAD source,
- reordering R function formals,
- editing R function bodies,
- adding gallery links to `.Rd` files,
- changing examples in a way that alters runtime or expected output,
- adding new method promises not already supported by the code.

## Current Inventory Signal

The largest flat argument surfaces are:

- `crs.Rd`: 59 argument items
- `clsd.Rd`: 40 argument items
- `krscvNOMAD.Rd`: 42 argument items
- `frscvNOMAD.Rd`: 35 argument items
- `krscv.Rd`: 25 argument items
- `crsiv.Rd`: 25 argument items
- `frscv.Rd`: 22 argument items
- `crsivderiv.Rd`: 19 argument items
- `crssigtest.Rd`: 18 argument items
- `snomadr.Rd`: 20 argument items

These are good candidates for argument subheaders because users currently have
to scan long historically ordered lists to find data inputs, spline controls,
cross-validation controls, NOMAD controls, output controls, and pass-through
options.

## Documentation Format Contract

Use real Rd subsections inside `\arguments{}`:

```r
\subsection{Model Inputs}{
  Core data and formula arguments.
}

\item{formula}{...}
\item{data}{...}
```

Rules:

- Keep top-level `\item{arg}{...}` entries top-level.
- Do not nest `\item{}` entries inside `\subsection{}` bodies.
- Do not create fake group headings as `\item{Group}{...}`.
- Do not create empty subsection bodies.
- Keep the visible order of important leading arguments intuitive for canonical
  use, especially `formula`, `data`, `xz`, `y`, `z`, `w`, and `x`.
- Organize by user mental model, not internal history.
- Preserve all existing aliases, usage signatures, examples, references, and
  value entries unless a documented mismatch is found.

## Proposed User-Mental-Model Groupings

### `crs.Rd`

Suggested argument groups:

- Model Inputs
- Basis, Spline, And Kernel Structure
- Cross-Validation And Search Controls
- NOMAD Controls
- Quantile, Weights, And Derivatives
- Returned State And Output Controls
- Warnings And Progress
- Pass-Through Arguments

Clarify that `crs.formula()` exposes the full search surface, while lower-level
helpers such as `frscv*()` and `krscv*()` document the exhaustive/NOMAD
subroutes in greater detail.

### `frscv.Rd` and `krscv.Rd`

Suggested groups:

- Data Inputs
- Basis And Spline Complexity
- Exhaustive Search Controls
- Kernel Or Factor Controls
- Quantile, Weights, And Numerical Guardrails
- Warnings And Progress

Make the distinction clear:

- `frscv()` searches factor/inclusion structure for categorical predictors.
- `krscv()` searches categorical-kernel bandwidths through `optim()`.

### `frscvNOMAD.Rd` and `krscvNOMAD.Rd`

Suggested groups:

- Data Inputs
- Basis And Spline Complexity
- NOMAD Search Controls
- Kernel Or Factor Controls
- Quantile, Weights, And Numerical Guardrails
- Warnings And Progress

Make the route-specific NOMAD default profiles easy to find, including the
fact that user-supplied `opts` entries take precedence.

### `snomadr.Rd`

Suggested groups:

- Objective, Dimension, And Evaluation Environment
- Starting Values, Bounds, And Variable Types
- Multi-Start Controls
- NOMAD Options And Information Queries
- Progress And Pass-Through Arguments

Preserve the NOMAD option reference link and compatibility-profile explanation,
but make it easier for a user to answer: "What do I pass for my problem?" and
"Where do NOMAD options go?"

### `crsiv.Rd` and `crsivderiv.Rd`

Suggested groups:

- Data Inputs
- Evaluation Inputs
- Regularization And Iteration Controls
- Derivatives And IV-Specific Controls
- CRS Pass-Through Search Controls
- NOMAD And Progress Controls
- Warnings And Output Controls

Clarify `...` as arguments passed through to `crs()`, and point users to
`crs.Rd` for the full spline/CV/NOMAD argument surface.

### `clsd.Rd`

Suggested groups:

- Data Inputs
- Density Basis Structure
- Constraint And Shape Controls
- Optimization Controls
- Evaluation, Plotting, And Output Controls
- Warnings, Messages, And Pass-Through Arguments

Because this page has a large argument surface and a specialized density
workflow, preserve exact semantics and avoid broad rewriting.

### Smaller Helper And Data Pages

Keep these mostly as-is unless they have clear documentation defects:

- helper pages with one or two arguments do not need subsection structure,
- data pages should remain lightweight and stable,
- do not add gallery links to `.Rd` files; the vignette is the intended pointer
  surface for gallery material.

## Implementation Protocol

1. Create a detached worktree under `/Users/jracine/Development/tmp`, not in the
   live `crs` tree the user may rebuild from.
2. In the detached worktree, run a baseline documentation and check snapshot:
   - `R CMD build`
   - `R CMD check --no-manual --ignore-vignettes`
   - `tools::checkRd()` across `man/*.Rd`
   - `R CMD Rd2txt` spot checks for the largest pages
3. Proof-of-concept one high-value page first, preferably `crs.Rd`.
4. Verify the proof-of-concept end to end:
   - `tools::checkRd("man/crs.Rd")`
   - `R CMD Rd2txt man/crs.Rd`
   - visual scan of rendered text for clean headers and no argument loss
   - package check smoke if necessary
5. Proceed in narrow clusters:
   - cluster 1: `crs.Rd`
   - cluster 2: `frscv*.Rd` and `krscv*.Rd`
   - cluster 3: `snomadr.Rd`
   - cluster 4: `crsiv*.Rd`
   - cluster 5: `clsd.Rd`
   - cluster 6: small helpers, package page, data pages only as needed
6. After each cluster:
   - compare argument names against function formals,
   - compare before/after `\usage{}` signatures,
   - run `tools::checkRd()` for touched pages,
   - render touched pages with `R CMD Rd2txt`,
   - record proof logs under `/Users/jracine/Development/tmp`.
7. After all clusters:
   - run full tests,
   - build and check a tarball,
   - install into a private library,
   - run representative installed smokes.
8. Merge/cherry-pick into the live `crs` repo only after green checks.
9. Commit in checkpoints so any problematic cluster can be reverted cleanly.

## Validation Gates

Documentation gates:

- `tools::checkRd()` passes for every touched page.
- `R CMD Rd2txt` output is readable for the major pages.
- `R CMD check --no-manual --ignore-vignettes` passes.
- No new undocumented argument, codoc, alias, or example warning appears.

API-preservation gates:

- All R function formal lists are unchanged in this documentation tranche.
- All R function bodies are unchanged in this documentation tranche.
- `\usage{}` entries remain synchronized with current formals.
- No examples are made slower or semantically different unless separately
  justified.
- `NAMESPACE`, `DESCRIPTION`, `R/`, `src/`, and `data/` remain untouched unless
  a documentation mismatch proves a tiny metadata edit is necessary.

Installed smoke gates:

- `library(crs)`
- `example(crs, run.dontrun = FALSE)`
- `example(snomadr, run.dontrun = FALSE)`
- selected small direct calls for:
  - `crs(..., cv="none")`
  - `crs(..., cv="nomad", max.bb.eval` small enough for smoke)
  - `frscv()` or `krscv()` with bounded small grids
  - `crsiv()` and `crsivderiv()` smoke using existing test fixtures if cheap
  - `clsd()` smoke using existing tests/examples if cheap

Regression gates:

- `testthat` suite passes.
- Demos are not edited in this tranche.
- Any failures are classified as harness/documentation/example failures before
  being attributed to package code.

## Critique From Highest Engineering Standards

The first naive plan would be to reorganize all `.Rd` files and function
formals together because that sounds like a complete cleanup. That is too risky.
R users may call functions positionally, and even an aesthetically superior
formal order can silently change meaning. Documentation order and R formal order
are separate risk axes and should not be mixed.

The second risk is over-broad rewriting. These pages contain statistical
semantics, NOMAD route details, and historical contract language. Rewriting
prose for style can accidentally weaken or change the estimator contract. The
campaign should move and clarify text first, not reinterpret it.

The third risk is trusting Rd syntax by inspection. The subsection pattern is
valid only if rendered and checked in the real package. Each cluster must run
`checkRd` and `Rd2txt`, and the final tarball check is the decisive proof.

The fourth risk is documentation drift in pass-through arguments. `crsiv()` and
`crsivderiv()` rely on `...` forwarding to `crs()`, while `crs()` delegates to
CV helpers. The docs should help users discover the deeper argument surface
without duplicating every downstream argument in every page.

The fifth risk is letting data/helper pages become noisier. Small pages with one
or two arguments do not benefit from grouped headings; applying the format
mechanically would reduce clarity.

## Reformulated Plan

The refined campaign is:

1. Make this a detached-worktree, documentation-only campaign.
2. Preserve R function formals exactly.
3. Prove the argument-subsection format on one major page before scaling.
4. Organize long argument lists by user mental model, with short, real
   subsection bodies and top-level argument items.
5. Clarify `...` and route delegation so users know when to consult `crs.Rd`,
   `snomadr.Rd`, or the CV helper pages.
6. Avoid broad prose rewrites; edit only where grouping, accuracy, or
   discoverability improves.
7. Validate each cluster with `checkRd`, `Rd2txt`, and argument/formal
   comparisons.
8. Run final tarball, private-library install, testthat, and installed smoke
   gates before live merge.
9. Do not pursue R formal reordering in this campaign; `.Rd` organization is
   sufficient for the current user-experience goal.

## Implementation Decision

Jeffrey confirmed that touching functions is unnecessary if documentation
reordering achieves the goal. Proceed with `.Rd` organization only: make the
manual pages more orderly, easier to scan, and more faithful to the current
code, without changing function formals, function bodies, defaults, examples,
or runtime behavior.
