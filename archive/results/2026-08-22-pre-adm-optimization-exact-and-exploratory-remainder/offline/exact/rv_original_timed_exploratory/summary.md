# RV-paper original timed solver: exploratory performance

Matrix: phi5 and phi6; durations 4, 8, and 16; epsilon in {1,2,4,8} when epsilon <= duration; samples 0, 1, and 2. The trace files are byte-identical to the archived `dataNew` inputs.

Each positive and negative check had a 5-second Z3 cap. Exact instance time is max(positive time, negative time). Solver/module import time is excluded. Four persistent worker processes were used.

The archived SMT assertions were not changed. The harness only replaced its Solver constructor with a timeout-aware wrapper that raises on Z3 `unknown`; this prevents the archived code from misrecording a timeout as an `unsat` proof.

## Overall

| Formula | Instances | Completed | Timed out | Completion | Median exact time (completed) | Maximum exact time (completed) |
|---|---:|---:|---:|---:|---:|---:|
| phi5 | 33 | 24 | 9 | 72.7% | 0.317 s | 5.289 s |
| phi6 | 33 | 21 | 12 | 63.6% | 0.256 s | 7.278 s |

Across both formulas, 45/66 instances completed and 21 timed out. There were 20 positive-check timeouts and 12 negative-check timeouts; 11 instances timed out in both polarities.

A reported wall time can exceed five seconds because the cap applies to Z3's `check()`, while Python constructs the explicit clock model before that call. This is most visible at d=16, eps=8.

## Comparison with random_grid_timed.py

| Outcome on same instance | Count |
|---|---:|
| Both completed | 13 |
| random_grid timed out; RV completed | 32 |
| random_grid completed; RV timed out | 1 |
| Both timed out | 20 |

The RV version completed 45/66 instances versus 14/66 for `random_grid_timed.py`. This is a timing comparison, not a semantic equivalence claim: the encodings differ, and the archived solver has known correctness defects.

For instances completed by both implementations:

| Formula | Common instances | random_grid median | RV median | Median per-instance random_grid/RV ratio |
|---|---:|---:|---:|---:|
| phi5 | 7 | 0.026 s | 0.551 s | 0.13x |
| phi6 | 6 | 0.038 s | 0.447 s | 0.15x |

## Completion matrix

Entries are completed/3 samples; a dash means epsilon exceeds duration.

| Formula | Duration | eps=1 | eps=2 | eps=4 | eps=8 |
|---|---:|---:|---:|---:|---:|
| phi5 | 4 | 3/3 | 3/3 | 3/3 | - |
| phi5 | 8 | 3/3 | 3/3 | 2/3 | 2/3 |
| phi5 | 16 | 3/3 | 2/3 | 0/3 | 0/3 |
| phi6 | 4 | 3/3 | 3/3 | 3/3 | - |
| phi6 | 8 | 3/3 | 3/3 | 2/3 | 1/3 |
| phi6 | 16 | 3/3 | 0/3 | 0/3 | 0/3 |

## Timed-out instances

### phi5

- d=8: eps=4: samples 0; eps=8: samples 1
- d=16: eps=2: samples 1; eps=4: samples 0,1,2; eps=8: samples 0,1,2

### phi6

- d=8: eps=4: samples 1; eps=8: samples 1,2
- d=16: eps=2: samples 0,1,2; eps=4: samples 0,1,2; eps=8: samples 0,1,2

Polarity details and raw elapsed times are in `timed_out_instances.csv`.

## Correctness caveat

These return values must not be used as a correctness oracle. Among other issues, the archived negative query uses a shared existential witness (effectively exists-v/forall-u instead of forall-u/exists-v), and its finite-horizon windows are boundary-sensitive. The experiment is retained only as a historical performance measurement.

Raw data: `results.csv`; machine-readable summary: `summary.json`; flat timeout list: `timed_out_instances.csv`.
