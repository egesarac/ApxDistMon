# Results

This directory separates the manuscript's reference measurements from new
benchmark output and derived figures.

- `reference/` contains the preserved result snapshot used for the paper.
- `online/` and `offline/` are created by new benchmark runs and may contain
  resumable SQLite checkpoints.
- `figures/` is created by the visualization scripts and contains derived
  PDFs.

The benchmark runners never write to `reference/` by default. Keep that
directory unchanged. To run an independent experiment, use the default runner
locations or provide a new `--output-dir`. To visualize the reference snapshot,
run the plotting commands in the main [benchmark README](../README.md).
