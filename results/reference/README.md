# Manuscript reference results

This directory is the preserved result snapshot supporting the manuscript.
It is input data for analysis and visualization, not an output directory for
new experiments.

## Layout

- `offline/approximate/` contains raw approximate-monitor timing and verdict
  rows for the random traces and both case studies.
- `offline/exact/` contains the matching raw exact-monitor timing and verdict
  rows.
- Exact random-trace results for timed formulas `phi5` and `phi6` were
  evaluated only at durations 4 and 8, as specified by the current experiment
  plan. Their duration-16 and duration-32 plot cells are intentionally blank;
  they are not missing or damaged data.
- `offline/timeouts.tsv` records the six random-trace and two
  mutual-separation exact-monitor timeouts used by the plots.
- `online/incremental/` and `online/naive/` contain raw online timing and
  verdict-sequence rows.
- `online/summary.csv` is the derived six-formula table consumed by the online
  plotter. The raw online directories also preserve the historical `ec` and
  `ed` files and durations 16 and 32; the default paper plot excludes those
  durations and uses only the six formulas in `summary.csv`.
- Generated PDFs are not stored in this reference snapshot. New plots are
  written outside this directory.

From the benchmark root (the directory containing `README.md` and
`CMakeLists.txt`), regenerate all figures with:

```bash
.venv/bin/python visualization/plot_online.py
.venv/bin/python visualization/plot_offline_random.py
.venv/bin/python visualization/plot_offline_ms.py
.venv/bin/python visualization/plot_offline_wt.py
```

The commands read this snapshot by default and write their output under
`results/figures/`. Use explicit `--input` or `--generated-input` arguments to
analyze a new run instead.

## Provenance and interpretation

This snapshot preserves materialized TXT, TSV, and CSV results. Historical
runner checkpoint databases, their WAL/SHM sidecars, and lock files are
intentionally omitted. Its rows were assembled through resumed and historical
runs; the retained artifacts do not establish that every row came from one
uninterrupted execution. Treat the timings as the manuscript's reference
measurements, not as bit-for-bit golden values or resumable runner state.

The six exact RG timeout jobs have no raw TXT rows. The offline RG plot uses
their nominal 600 s limits and six explicit verdicts independently established
by the [archived continuous real-clock EDM
run](../../archive/results/2026-08-22-pre-adm-optimization-selective-reset/).
The two exact mutual-separation timeouts likewise have no raw rows and are
drawn as crosses at 360 s. These timeout substitutions are part of the
analysis provenance and must not be inferred from the TXT files alone.

The online raw files and `summary.csv` are a historical superset of the
current runner plan: they retain durations 16 and 32 and legacy formula IDs
`ec` and `ed`. The default paper plot filters to durations 64--512 and the six
formula IDs in the current plan. See [`../../REPRODUCING.md`](../../REPRODUCING.md)
for current-grid counts and acceptance criteria.
