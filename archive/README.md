# Code archive

This directory contains old work from after the original STTT submission. The
original RV and STTT artifacts are preserved in the original monorepo's pinned
[top-level `archive/`
tree](https://github.com/egesarac/distributed-monitoring/tree/f68294625456612c3a1ece00439ebc428a67eb01/archive);
that companion tree is not part of a flattened standalone checkout. Files here
are for reference and are not current benchmark results.

## Contents

- `analysis/` contains old notes and diagnostic programs. Superseded drafts
  are clearly named under `analysis/editor-snapshots/`. Some paths may be stale.
- `solvers/2026-08-22-timed-edm-real-clock-and-discrete-grid-experiments/`
  contains experimental timed EDM solvers. The continuous solver uses
  real-valued clocks. The discrete-grid solver has different semantics. None
  is an active backend.
- `source/retired-monitoring-headers/` contains old headers that active code
  does not include.
- `figures/offline-random/` contains the two cross-snapshot August figures.
  Other figures stay with their result snapshot. August 22 uses the older ADM.
  August 23 uses the optimized ADM, especially for untimed formulas. Both use
  the same optimized real-clock EDM data.

## Result snapshots

- `results/2026-08-22-pre-adm-optimization-selective-reset/` contains
  older ADM results. Timed EDM for `phi5` and `phi6` already used the optimized
  real-clock solver.
- `results/2026-08-22-pre-adm-optimization-exact-and-exploratory-remainder/`
  contains untimed EDM results and incomplete timed experiments.
- `results/2026-08-23-post-adm-optimization-generated-results/`
  contains optimized ADM results. Its timed EDM data was copied from August 22
  for the historical figure; it is not a new timed run.
- `results/2026-08-24-case-study-rerun-history/` mixes complete,
  preliminary, and interrupted runs.
- `results/2026-08-26-invalid-timed-grid-scale-2eps/` contains timed
  results that are invalid when `d != eps`. Only rows with `d = eps` are valid.
  The active solver uses `2*d`; see the
  [counterexample](analysis/scale_2eps_plus_1_counterexample.md) and
  [safety note](analysis/timed_grid_scale_2d_safety.md).

Historical runner checkpoint databases, their WAL/SHM sidecars, and lock files
are intentionally omitted. The materialized result files and derived figures
are the archival record; these snapshots are not resumable runner outputs. Do
not merge them into `results/`.
