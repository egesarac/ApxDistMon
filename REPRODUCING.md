# Reproducing the paper experiments

All commands and relative paths below start at the **benchmark root**, the
repository root containing `README.md` and `CMakeLists.txt`. The default
runners write a fresh, resumable run to `results/offline/` and
`results/online/`. The shipped measurements under `results/reference/` are
read-only inputs to the plotting scripts; they are not runner checkpoints.

The supported and reference platform, dependency installation, build, test,
runner, and plotting commands are in [`README.md`](README.md). Input and
reference-result provenance are documented in [`data/README.md`](data/README.md)
and [`results/reference/README.md`](results/reference/README.md).

## Paper-to-code map

The paper and runners use the following formula identifiers:

| Paper | Runner ID | Formula |
| --- | --- | --- |
| `phi_1` | `phi1` | `G(p AND q)` |
| `phi_2` | `phi2` | `G(p OR q)` |
| `phi_3` | `phi3` | `p U q` |
| `phi_4` | `phi4` | `G(p => F q)` |
| `phi_5` | `phi5` | `G(p => F_[0,1) q)` |
| `phi_6` | `phi6` | `G(p => F_[0,2) q)` |
| `psi_1` | `ac` | past-always `(p AND q)` |
| `psi_2` | `ad` | past-always `(p OR q)` |
| `psi_3` | `since` | `p` since `q` |
| `psi_4` | `response` | past response |
| `psi_5` | `response01` | past response with `[0,1)` |
| `psi_6` | `response02` | past response with `[0,2)` |

The execution and artifact map is:

| Paper result | Runner stage and implementation | Fresh raw result |
| --- | --- | --- |
| Offline RG, `phi1`-`phi6` | `approximate-random`, CMake target `sttt_offline_approximate_random`; `exact-random`, `benchmarks/offline/exact/run.py random` | `results/offline/{approximate,exact}/random/phiN.txt`, `results/offline/timeouts.tsv` |
| Online RG, `psi1`-`psi6` | `incremental/ID` and `naive/ID`, CMake targets `sttt_online_MODE_ID` | `results/online/incremental/onlineIncr_ID.txt`, `results/online/naive/onlineNaive_ID.txt`, derived `results/online/summary.csv` |
| Water tank, `phi_WT = G(sum(x_i) > 10)` | `approximate-water-tanks`, target `sttt_offline_approximate_water_tanks`; `exact-water-tanks`, `run.py water-tanks` | `results/offline/approximate/case_studies/water_tanks*.txt`, `results/offline/exact/case_studies/water_tanks.txt` |
| Swarm of drones, `phi_SD` (mutual separation) | `approximate-mutual-separation`, target `sttt_offline_approximate_mutual_separation`; `exact-mutual-separation`, `run.py mutual-separation` | `results/offline/approximate/case_studies/mutual_separation*.txt`, `results/offline/exact/case_studies/mutual_separation.txt` |

For case-study approximate results, no suffix means ADM, `_adm_f` means
ADM-F, `_adm_fr` means ADM-FR, `_adm_c` means ADM-C, and `_adm_cr` means
ADM-CR. The coarse variants are not applicable to mutual separation.

For the shipped snapshot, insert `reference/` after `results/` in the raw
paths above. The paper-figure map is:

| Manuscript label and asset | Plot command | Generated output |
| --- | --- | --- |
| `fig:rgresults`, `speedup_all.pdf` | `.venv/bin/python visualization/plot_offline_random.py` | `results/figures/speedup_offline_random.pdf` |
| `fig:rgresultsOnline`, `speedup_online.pdf` | `.venv/bin/python visualization/plot_online.py` | `results/figures/speedup_online.pdf` |
| `fig:wtresults`, `speedup_offline_wt.pdf` | `.venv/bin/python visualization/plot_offline_wt.py` | `results/figures/speedup_offline_wt.pdf` |
| `fig:sdresults`, `speedup_offline_ms.pdf` | `.venv/bin/python visualization/plot_offline_ms.py` | `results/figures/speedup_offline_ms.pdf` |

The plotters write only to `results/figures/` by default.

## Experiment contract

`scripts/run_offline.py --dry-run` and `scripts/run_online.py --dry-run` are
the authoritative executable plans. Their paper-default matrices are:

| Stage | Parameter grid | Repetitions | Timeout per solver check | Jobs |
| --- | --- | ---: | ---: | ---: |
| Approximate RG | `phi1..phi6`; `d={4,8,16,32}` s; `eps={1,2,4,8}` s with `eps<=d`; samples `0..99` | 100 | none | 9,000 |
| Exact RG | `phi1..phi4` on the same grid; `phi5,phi6` only at `d={4,8}`; samples `0..99` | 3 | 600 s | 7,400 |
| Approximate drones | `n={2,3,4}`; epsilon multiplier `{1,2,3,4,5}`; methods ADM, ADM-F, ADM-FR | 10 | none | 45 |
| Approximate water tanks | `n={2,3,4}`; epsilon multiplier `{1,2,4,8}`; methods ADM, ADM-F, ADM-FR, ADM-C, ADM-CR | 10 | none | 60 |
| Exact drones | `n={2,3,4}`; epsilon multiplier `{1,2,3,4,5}` | 10 | 360 s | 15 |
| Exact water tanks | `n={2,3,4}`; epsilon multiplier `{1,2,4,8}` | 10 | 120 s | 12 |
| Online RG | modes `{incremental,naive}`; six online IDs; `d={64,128,256,512}` s; `eps={2,4,8}` s; `Delta={4,8,16,32}` s with `Delta>eps` and `Delta<=d`; samples `0..99` | 10 | none | 43,200 |

The offline total is 16,532 jobs. For each RG sample `c`, the two signal
components are files `d_c.txt` and `d_(c+100).txt`. A case-study epsilon
multiplier is multiplied by the 0.05 s sampling period. The timed exact RG
restriction produces 74 paired `(formula,d,eps)` cells and 7,400 paired trace
instances; the additional approximate `phi5`/`phi6` results at `d=16,32` are
not plotted. Each online mode/formula stage has 3,600 jobs, for 43,200 total.

## Result contract

Raw TXT rows are whitespace-separated and have no header:

```text
approximate RG:
d_s eps_s threshold sample segments mean_seconds verdict

exact RG:
d_s eps_s - sample - positive_mean_seconds negative_mean_seconds positive_check negative_check verdict

approximate case study:
d_s n epsilon_multiplier segments mean_seconds false_possible true_possible verdict

exact case study:
d_s n epsilon_multiplier positive_mean_seconds positive_check negative_mean_seconds negative_check verdict

online:
d_s eps_s Delta_s sample mean_total_seconds verdict_sequence
```

`threshold` is the literal `0`, and the two exact-RG placeholder fields are
the literal `-`. All timing columns are arithmetic means over the repetitions
in the experiment table. Each exact-RG check time includes preprocessing. The
paper and case/offline plots report exact-monitor time as the larger of the
positive and negative mean check times, modeling the two checks in parallel.

Verdicts use `0 = false`, `1 = true`, and `2 = inconclusive (?)`. The exact
check flags are the Boolean results of the positive and negative proof checks;
the exact verdict is 1 for positive-only, 0 for negative-only, and 2 otherwise.
An online verdict sequence contains one digit per portion, in chronological
order, and has length `ceil(d/Delta)`.

The derived online CSV has this header:

```text
formula,d,eps,Delta,speedup,incr wall time,time per portion
```

For every complete 100-sample cell, `speedup` is the mean of the 100
per-sample naive/incremental ratios, `incr wall time` is the mean incremental
time in seconds, and `time per portion` divides that mean by `d/Delta`.
Incomplete cells are omitted. Materialization rejects an incremental/naive
verdict-sequence disagreement.

For each complete offline RG cell, the plotter averages the 100 trace times,
computes raw speedup as `sum(EDM) / sum(ADM)`, and reports the fraction where
ADM is inconclusive but EDM is conclusive. The combined-monitor cost is ADM
plus EDM only on traces where ADM is inconclusive. The plotter rejects any
conclusive ADM verdict that differs from EDM.

If any exact solver check reaches its limit, Z3 `unknown` becomes process exit
124. The runner records a terminal `timeout`, writes no raw result row for the
job, and materializes `timeouts.tsv` with columns
`stage, instance, timeout_seconds, reason`. A full offline run may therefore
exit successfully with timeouts; inspect `--status` and the counts below.

In the reference plots, each of the six RG timeouts is assigned 600 s and an
explicitly verified verdict recorded in `plot_offline_random.py`. The two
drone timeouts remain missing exact rows and are drawn as crosses at 360 s.
The water-tank and online experiments have no reference timeouts.

## Full-run checks and cost

Before a long run, verify the plan:

```bash
.venv/bin/python scripts/run_offline.py --dry-run
.venv/bin/python scripts/run_online.py --dry-run
```

The paper-default totals are 16,532 offline jobs and 43,200 online jobs. The
preserved snapshot contains 16,524 completed offline rows plus eight timeout
records: six exact `phi6` RG jobs and two exact drone jobs (`n=4`, epsilon
0.20 and 0.25 s). Together they cover all paper-default offline jobs. The
expected fresh raw counts for that outcome are:

- approximate RG: 1,500 rows in each of six files;
- exact RG: 1,500 rows for each of `phi1`-`phi4`, 700 for `phi5`, and 694 plus
  six timeout records for `phi6`;
- approximate drones: 15 rows per method file; approximate water tanks: 12
  rows per method file;
- exact drones: 13 rows plus two timeout records; exact water tanks: 12 rows;
- online: 3,600 rows in each of 12 current raw files and 216 data rows in
  `summary.csv`.

Timeout boundaries and timings are hardware-sensitive. For a fresh run on a
different machine, the deterministic acceptance criteria are that every
planned job is terminal, no job has status `error`, online paired verdicts
agree, no conclusive offline ADM verdict differs from EDM, the raw counts plus
timeout counts cover the plan, and all four plotting commands exit successfully
and create nonempty PDFs. Explain any timeout-count deviation rather than
comparing timings bit-for-bit.

The runners execute benchmark jobs sequentially. Summing the stored mean
times times the paper repetitions, and adding one nominal limit per recorded
timeout, gives an approximate budget of 69.1 hours offline and 3.69 hours
online on the paper machine. These are not end-to-end measurements of one
uninterrupted run: they exclude build, process, and checkpoint overhead, and
the offline estimate only approximates preprocessing and timeout work.
