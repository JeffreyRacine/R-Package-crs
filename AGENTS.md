# AGENTS.md

Use `/Users/jracine/Development/AGENTS.md` as the canonical workflow and benchmarking policy.
Hard default: enforce `Sweep Safety Protocol (Default Across Scope Repos)` from the canonical AGENTS file.

Repo-specific note:
- Keep ridge/CV performance changes numerically stable across fit, predict, and CV pathways.
- Native-code release-hardening reminder: for changes touching `src/`, NOMAD
  interface code, registered native interfaces, or `.C`/`.Call` payload
  lifetimes, run the shared release gate with `RUN_RCHK=1` when container
  infrastructure is available:
  `cd /Users/jracine/Development && RUN_RCHK=1 ./release_protocol/run_crs_release_gate.sh`.
  `RUN_RCHK=auto` is acceptable for ordinary rehearsal only if the resulting
  summary records either PASS or a precise SKIP reason.

Documentation note:
- For long function `.Rd` argument lists, prefer real `\subsection{...}{...}`
  headers inside `\arguments{}` followed by top-level `\item{arg}{...}`
  entries. Organize arguments by user mental model while preserving the actual
  R function signatures and behavior.
