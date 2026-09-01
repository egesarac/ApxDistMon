# random_grid_timed exploratory performance

Matrix: phi5 and phi6; durations 4, 8, and 16; epsilon in {1,2,4,8} when epsilon <= duration; samples 0, 1, and 2.

Each positive and negative check had a 5-second Z3 cap. Exact instance time is max(positive time, negative time). Solver/module import time is excluded. Four persistent worker processes were used.

## Overall

| Formula | Instances | Completed | Timed out | Completion | Median exact time (completed) | Maximum exact time (completed) |
|---|---:|---:|---:|---:|---:|---:|
| phi5 | 33 | 7 | 26 | 21.2% | 0.026 s | 0.074 s |
| phi6 | 33 | 7 | 26 | 21.2% | 0.062 s | 1.256 s |

All 66 positive checks completed. Only 14 of 66 negative checks completed; every instance timeout came from the negative polarity.

Completed verdicts: phi5 had seven false verdicts. Phi6 had five false, one true, and one inconclusive verdict.

Four representative timeouts were retried with 30-second caps and all four timed out. The smallest (`phi5`, d=4, epsilon=1, sample=1, negative) was also retried alone and still timed out after 30 seconds.

## Completion matrix

Entries are completed/3 samples; a dash means that epsilon exceeds duration.

| Formula | Duration | eps=1 | eps=2 | eps=4 | eps=8 |
|---|---:|---:|---:|---:|---:|
| phi5 | 4 | 2/3 | 0/3 | 0/3 | - |
| phi5 | 8 | 1/3 | 0/3 | 1/3 | 0/3 |
| phi5 | 16 | 2/3 | 1/3 | 0/3 | 0/3 |
| phi6 | 4 | 2/3 | 1/3 | 0/3 | - |
| phi6 | 8 | 0/3 | 1/3 | 0/3 | 0/3 |
| phi6 | 16 | 2/3 | 1/3 | 0/3 | 0/3 |

## Timed-out instances

### phi5

- d=4: eps=1: samples 1; eps=2: samples 0,1,2; eps=4: samples 0,1,2
- d=8: eps=1: samples 0,1; eps=2: samples 0,1,2; eps=4: samples 0,1; eps=8: samples 0,1,2
- d=16: eps=1: samples 1; eps=2: samples 0,1; eps=4: samples 0,1,2; eps=8: samples 0,1,2

### phi6

- d=4: eps=1: samples 2; eps=2: samples 0,2; eps=4: samples 0,1,2
- d=8: eps=1: samples 0,1,2; eps=2: samples 0,1; eps=4: samples 0,1,2; eps=8: samples 0,1,2
- d=16: eps=1: samples 1; eps=2: samples 0,1; eps=4: samples 0,1,2; eps=8: samples 0,1,2

Every listed timeout occurred in the negative check.

Raw data: `results.csv`; machine-readable summary: `summary.json`; flat timeout list: `timed_out_instances.csv`; 30-second retries: `retries_30s.csv`.
