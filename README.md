# Benchmark code

This directory contains the maintained monitors, benchmarks, tests, input
data, and figure scripts. In this documentation, the **benchmark root** is the
repository root, which contains this `README.md` and `CMakeLists.txt`. Run the
commands below from the benchmark root.

You need CMake 3.16 or newer and Python 3.12. The maintained
benchmark runners require a POSIX system providing `fcntl.flock`, and the C++
implementation requires GCC or Clang on a target supporting `__int128`.
Windows/MSVC is not supported. The paper's reference environment was Ubuntu
24.04 on x86-64, g++ 13.3.0, and Python 3.12.3.

## Set up Python

```bash
python3.12 -m venv .venv
.venv/bin/pip install -r requirements-smt.txt -r requirements-plotting.txt
```

## Build and test

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure
.venv/bin/python -m unittest discover -s tests -p 'test_*.py'
```

## Run benchmarks

```bash
.venv/bin/python scripts/run_online.py
.venv/bin/python scripts/run_offline.py
```

New results and checkpoints go to `results/online/` and `results/offline/`.
These directories are separate from the shipped manuscript results in
`results/reference/`. Interrupted new runs resume from their checkpoints. Use
`--list` to see available experiments, `--dry-run` to inspect the selected
plan, and `--help` to see filters. Do not use `--restart` with a reference
result path.

## Results layout

- `results/reference/` is the read-only snapshot supporting the manuscript.
- `results/online/` and `results/offline/` are created by new benchmark runs.
- `results/figures/` contains derived plots and can be regenerated.

See [`results/README.md`](results/README.md) for the detailed layout.
See [`REPRODUCING.md`](REPRODUCING.md) for the paper-to-code map, full
experiment and result contracts, expected costs, and completion criteria.

## Make figures

The following commands read `results/reference/` by default and write PDFs to
`results/figures/`:

```bash
.venv/bin/python visualization/plot_online.py
.venv/bin/python visualization/plot_offline_random.py
.venv/bin/python visualization/plot_offline_ms.py
.venv/bin/python visualization/plot_offline_wt.py
```

To plot a new run instead, pass `--input results/online/summary.csv` to the
online plotter or `--generated-input results/offline` to an offline plotter.
Use `--out` to choose a different figure path.

## License

The contents of this directory are available under the [MIT License](LICENSE).
