# Vendored NOMAD Patch Manifest

Last updated: 2026-07-16

## Upstream baseline

`crs` vendors NOMAD 4.5.0 under `src/nomad4_src`. The import is recorded by
commit `e872726` (`Upgrade crs NOMAD bridge to embedded NOMAD 4.5.0`,
2026-02-23), and the embedded version remains defined in
`src/nomad4_src/src/nomad_version.hpp`.

This directory is not an unmodified upstream copy. Future NOMAD upgrades must
review and deliberately preserve, replace, or retire the local adaptations
below rather than overwriting the tree wholesale.

## Local adaptation groups

1. **R-native lifecycle and C interface**
   - exception-path cleanup and reusable global NOMAD state (`ed3a0b3`);
   - cache/evaluation accounting exposed through the C interface (`8720cd2`);
   - versioned observed-solve and cooperative-poll support, bounded callback
     output handling, and explicit interrupt ownership (`129dcac`);
   - signal-safe SIGINT capture, orderly checkpoint consumption, handler
     restoration, interrupt precedence, and retained R-callback diagnostics
     (`528c9e8`).
2. **R console and failure integration**
   - NOMAD and SGTELIB diagnostics are routed through R console functions
     instead of raw process stdout/stderr or process-level exits. The March 23,
     2026 commit series beginning with `cd6200a` and ending with `307bac5`
     records these changes.
3. **Random-number reproducibility**
   - SGTELIB surrogate randomness uses NOMAD's controlled RNG rather than the
     process-global system `rand()` (`4572034`).
4. **Portability and build contracts**
   - compiler-warning repairs and floating-point qualification (`9e4cdb0`,
     `6f055d8`, `c978bff`);
   - the shipped build intentionally preserves the package's sealed OpenMP
     contract (`8c5825a`).

The authoritative file-level record is Git history:

```sh
git log -- src/nomad4_src
```

## Upgrade checklist

For any future upstream NOMAD refresh:

1. record the exact upstream version/commit and license state;
2. compare upstream against the current vendored tree before replacing files;
3. classify every adaptation group above as retained, superseded, or retired;
4. preserve registered callable names, API versions, enum values, struct
   layouts, callback ownership, and successful-run numerical/accounting
   behavior unless an explicitly approved API tranche changes them;
5. revalidate exception cleanup, consecutive solves, real SIGINT and repeated
   SIGINT, R-handler restoration, observer stop behavior, callback diagnostics,
   RNG reproducibility, console routing, and cache/evaluation accounting;
6. run focused downstream installed-tarball proof in `np` and `npRmpi`,
   including MPI external-interrupt ownership, before broader release gates;
7. run compiler-warning gates and Docker `rchk` on the resulting `crs`
   tarball.

Do not resolve an upstream merge by silently dropping a local adaptation. If a
local patch is no longer needed, document the upstream replacement and the
validation that permits its retirement.
