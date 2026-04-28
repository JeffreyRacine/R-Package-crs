# AGENTS.md

Use `/Users/jracine/Development/AGENTS.md` as the canonical workflow and benchmarking policy.
Hard default: enforce `Sweep Safety Protocol (Default Across Scope Repos)` from the canonical AGENTS file.

Repo-specific note:
- Keep ridge/CV performance changes numerically stable across fit, predict, and CV pathways.

Documentation note:
- For long function `.Rd` argument lists, prefer real `\subsection{...}{...}`
  headers inside `\arguments{}` followed by top-level `\item{arg}{...}`
  entries. Organize arguments by user mental model while preserving the actual
  R function signatures and behavior.
