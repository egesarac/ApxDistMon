#!/usr/bin/env python3
"""Generate the offline random-trace benchmark heatmaps."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LogNorm


BENCHMARK_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_GENERATED_INPUT = BENCHMARK_ROOT / "results/reference/offline"
DEFAULT_OUTPUT = BENCHMARK_ROOT / "results/figures/speedup_offline_random.pdf"
INPUT_FILES = (
    "gand.csv",
    "gor.csv",
    "until.csv",
    "gimpf.csv",
    "gimpft1.csv",
    "gimpft2.csv",
)
FORMULAS = tuple(f"phi{index}" for index in range(1, 7))
DURATIONS = (4, 8, 16, 32)
EPSILON_VALUES = (1, 2, 4, 8)
SAMPLES_PER_CELL = 100
# The timeout report has no verdict column. These verdicts were independently
# established by the archived continuous real-clock EDM run in
# archive/results/2026-08-22-pre-adm-optimization-selective-reset.
# Keeping the list explicit prevents a future timeout from being silently
# classified as having the same verdict as ADM.
VERIFIED_EXACT_TIMEOUT_VERDICTS = {
    ("phi6", 8, 1, 12): 0,
    ("phi6", 8, 1, 35): 2,
    ("phi6", 8, 1, 48): 0,
    ("phi6", 8, 1, 79): 0,
    ("phi6", 8, 1, 97): 2,
    ("phi6", 8, 2, 48): 2,
}
DTYPE = {
    "d": int,
    "eps": int,
    "ADM": float,
    "EDM": float,
    "FP": float,
    "speedup": float,
    "label": str,
    "combined": float,
}


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot heatmaps from the offline random benchmark CSV files."
    )
    source = parser.add_mutually_exclusive_group()
    source.add_argument(
        "--input-dir",
        type=Path,
        help="directory containing archived input CSV files",
    )
    source.add_argument(
        "--generated-input",
        type=Path,
        help=(
            "result root containing approximate/random "
            f"and exact/random (default: {DEFAULT_GENERATED_INPUT})"
        ),
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"output figure path (default: {DEFAULT_OUTPUT})",
    )
    return parser.parse_args(argv)


def load_archived_results(input_dir: Path) -> pd.DataFrame:
    input_paths = [input_dir.expanduser().resolve() / name for name in INPUT_FILES]
    missing = [path for path in input_paths if not path.is_file()]
    if missing:
        missing_names = ", ".join(str(path) for path in missing)
        raise FileNotFoundError(f"missing offline result CSV files: {missing_names}")

    frames: list[pd.DataFrame] = []
    for formula, input_path in zip(FORMULAS, input_paths):
        frame = pd.read_csv(input_path, dtype=DTYPE)
        frame.insert(0, "formula", formula)
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def _result_rows(path: Path, expected_columns: int) -> list[list[str]]:
    if not path.is_file():
        raise FileNotFoundError(f"missing generated result: {path}")
    rows = [
        line.split()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    invalid = [
        index
        for index, row in enumerate(rows, start=1)
        if len(row) != expected_columns
    ]
    if invalid:
        raise ValueError(
            f"{path} has rows with the wrong number of fields: {invalid[:5]}"
        )
    return rows


def _nonnegative_time(raw_value: str, path: Path) -> float:
    value = float(raw_value)
    if not math.isfinite(value) or value < 0:
        raise ValueError(f"{path} contains invalid time {raw_value!r}")
    return value


def _verdict(raw_value: str, path: Path) -> int:
    value = int(raw_value)
    if value not in (0, 1, 2):
        raise ValueError(f"{path} contains invalid verdict {raw_value!r}")
    return value


def _format_combined(value: float) -> str:
    if not math.isfinite(value):
        return ""
    if value < 1:
        return "<1"
    return f"{value:.2g}"


def _verified_exact_timeout_results(
    input_path: Path,
) -> dict[tuple[str, int, int, int], tuple[float, int]]:
    report_path = input_path / "timeouts.tsv"
    if not report_path.is_file():
        return {}

    with report_path.open(encoding="utf-8", newline="") as report:
        reader = csv.DictReader(report, delimiter="\t")
        required_fields = {"stage", "instance", "timeout_seconds"}
        if not required_fields.issubset(reader.fieldnames or ()):
            raise ValueError(f"{report_path} has an invalid header")

        results: dict[tuple[str, int, int, int], tuple[float, int]] = {}
        for row in reader:
            if row["stage"] != "exact-random":
                continue
            instance = row["instance"].split("|")
            if len(instance) < 5 or instance[0] != "exact-random":
                raise ValueError(
                    f"{report_path} has an invalid exact-random instance: "
                    f"{row['instance']!r}"
                )
            key = (
                instance[1],
                int(instance[2]),
                int(instance[3]),
                int(instance[4]),
            )
            verdict = VERIFIED_EXACT_TIMEOUT_VERDICTS.get(key)
            if verdict is None:
                continue
            if key in results:
                raise ValueError(f"duplicate exact timeout for {key}")
            results[key] = (
                _nonnegative_time(row["timeout_seconds"], report_path),
                verdict,
            )
    return results


def load_generated_results(input_path: Path) -> pd.DataFrame:
    input_path = input_path.expanduser().resolve()
    timeout_results = _verified_exact_timeout_results(input_path)
    result_rows: list[dict[str, object]] = []

    for formula in FORMULAS:
        approximate_path = input_path / "approximate/random" / f"{formula}.txt"
        exact_path = input_path / "exact/random" / f"{formula}.txt"
        approximate: dict[tuple[int, int, int], tuple[float, int]] = {}
        exact: dict[tuple[int, int, int], tuple[float, int]] = {}
        imputed_timeouts: set[tuple[int, int, int]] = set()

        for fields in _result_rows(approximate_path, 7):
            key = (int(fields[0]), int(fields[1]), int(fields[3]))
            if key in approximate:
                raise ValueError(f"duplicate approximate result for {formula} {key}")
            approximate[key] = (
                _nonnegative_time(fields[5], approximate_path),
                _verdict(fields[6], approximate_path),
            )

        for fields in _result_rows(exact_path, 10):
            key = (int(fields[0]), int(fields[1]), int(fields[3]))
            if key in exact:
                raise ValueError(f"duplicate exact result for {formula} {key}")
            exact[key] = (
                max(
                    _nonnegative_time(fields[5], exact_path),
                    _nonnegative_time(fields[6], exact_path),
                ),
                _verdict(fields[9], exact_path),
            )

        for (
            timeout_formula,
            duration,
            epsilon,
            sample,
        ), result in timeout_results.items():
            if timeout_formula != formula:
                continue
            key = (duration, epsilon, sample)
            if key not in exact:
                exact[key] = result
                imputed_timeouts.add(key)

        for duration in DURATIONS:
            for epsilon in EPSILON_VALUES:
                valid = epsilon <= duration
                keys = [
                    (duration, epsilon, sample)
                    for sample in range(SAMPLES_PER_CELL)
                ]
                paired = [key for key in keys if key in approximate and key in exact]
                complete = valid and len(paired) == SAMPLES_PER_CELL
                row: dict[str, object] = {
                    "formula": formula,
                    "d": duration,
                    "eps": epsilon,
                    "paired_samples": len(paired),
                    "imputed_timeouts": sum(
                        key in imputed_timeouts for key in keys
                    ),
                    "complete": complete,
                    "ADM": float("nan"),
                    "EDM": float("nan"),
                    "FP": float("nan"),
                    "speedup": float("nan"),
                    "combined": float("nan"),
                }

                if complete:
                    approximate_times = [approximate[key][0] for key in keys]
                    exact_times = [exact[key][0] for key in keys]
                    approximate_verdicts = [approximate[key][1] for key in keys]
                    exact_verdicts = [exact[key][1] for key in keys]
                    unsound = [
                        key
                        for key, approximate_verdict, exact_verdict in zip(
                            keys, approximate_verdicts, exact_verdicts
                        )
                        if approximate_verdict != 2
                        and approximate_verdict != exact_verdict
                    ]
                    if unsound:
                        raise ValueError(
                            f"{formula} has inconsistent conclusive verdicts: "
                            f"{unsound[:5]}"
                        )

                    approximate_total = sum(approximate_times)
                    exact_total = sum(exact_times)
                    combined_total = sum(
                        approximate_time
                        + (exact_time if approximate_verdict == 2 else 0.0)
                        for approximate_time, exact_time, approximate_verdict in zip(
                            approximate_times,
                            exact_times,
                            approximate_verdicts,
                        )
                    )
                    row.update(
                        {
                            "ADM": approximate_total / SAMPLES_PER_CELL,
                            "EDM": exact_total / SAMPLES_PER_CELL,
                            "FP": 100.0
                            * sum(
                                approximate_verdict == 2 and exact_verdict != 2
                                for approximate_verdict, exact_verdict in zip(
                                    approximate_verdicts, exact_verdicts
                                )
                            )
                            / SAMPLES_PER_CELL,
                            "speedup": exact_total / approximate_total,
                            "combined": exact_total / combined_total,
                        }
                    )
                result_rows.append(row)

    return pd.DataFrame(result_rows)


def plot_results(data: pd.DataFrame, output_path: Path) -> None:
    speedup_values = data["speedup"].to_numpy(dtype=float)
    speedup_values = speedup_values[
        np.isfinite(speedup_values) & (speedup_values > 0)
    ]
    combined_values = data["combined"].to_numpy(dtype=float)
    combined_values = combined_values[
        np.isfinite(combined_values) & (combined_values > 0)
    ]
    if speedup_values.size == 0 or combined_values.size == 0:
        raise ValueError("offline results contain no positive speedup values")
    speedup_norm = LogNorm(
        vmin=float(speedup_values.min()),
        vmax=float(speedup_values.max()),
    )
    combined_norm = LogNorm(
        vmin=float(combined_values.min()),
        vmax=float(combined_values.max()),
    )

    figure, axes = plt.subplots(len(FORMULAS), 3)
    figure.set_size_inches([12, 12.5])

    for row_index, formula in enumerate(FORMULAS):
        formula_data = data[data["formula"] == formula]
        false_positives = formula_data.pivot_table(
            "FP", "eps", "d", fill_value=float("nan"), dropna=False
        ).reindex(index=EPSILON_VALUES, columns=DURATIONS)
        combined = formula_data.pivot_table(
            "combined", "eps", "d", fill_value=float("nan"), dropna=False
        ).reindex(index=EPSILON_VALUES, columns=DURATIONS)
        combined_annotations = combined.apply(
            lambda column: column.map(_format_combined)
        )
        speedup = formula_data.pivot_table(
            "speedup", "eps", "d", fill_value=float("nan"), dropna=False
        ).reindex(index=EPSILON_VALUES, columns=DURATIONS)

        sns.heatmap(
            speedup,
            cmap="GnBu",
            annot=speedup,
            fmt=".2g",
            linewidth=0.5,
            ax=axes[row_index][0],
            norm=speedup_norm,
        )
        axes[row_index][0].set(xlabel="duration", ylabel="ε")
        axes[row_index][0].set_yticks(
            [0.5, 1.5, 2.5, 3.5], labels=list(speedup.index.values)
        )

        sns.heatmap(
            false_positives,
            cmap="GnBu",
            annot=false_positives,
            fmt=".0f",
            linewidth=0.5,
            ax=axes[row_index][1],
            vmin=0,
            vmax=100,
        )
        axes[row_index][1].set(xlabel="duration", ylabel="")
        axes[row_index][1].set_yticks(
            [0.5, 1.5, 2.5, 3.5], labels=list(false_positives.index.values)
        )

        sns.heatmap(
            combined,
            cmap="GnBu",
            annot=combined_annotations,
            fmt="",
            linewidth=0.5,
            ax=axes[row_index][2],
            norm=combined_norm,
        )
        axes[row_index][2].set(xlabel="duration", ylabel="")
        axes[row_index][2].set_yticks(
            [0.5, 1.5, 2.5, 3.5], labels=list(false_positives.index.values)
        )

        if row_index == 0:
            axes[row_index][0].set_xticks(
                [0.5, 1.5, 2.5, 3.5], labels=speedup.columns
            )
            axes[row_index][0].set_title("Speed-up (times) ADM vs. EDM")
            axes[row_index][1].set_xticks(
                [0.5, 1.5, 2.5, 3.5], labels=false_positives.columns
            )
            axes[row_index][1].set_title("Extra Inconclusive Rate (%)")
            axes[row_index][2].set_xticks(
                [0.5, 1.5, 2.5, 3.5], labels=false_positives.columns
            )
            axes[row_index][2].set_title("Speed-up (times) Combined vs. EDM")

    row_labels = [rf"$\varphi_{number}$" for number in range(1, 7)]
    for axis, row_label in zip(axes[:, 0], row_labels):
        axis.annotate(
            row_label,
            xy=(0, 0.5),
            xytext=(-axis.yaxis.labelpad - 5, 0),
            xycoords=axis.yaxis.label,
            textcoords="offset points",
            size="large",
            ha="right",
            va="center",
            rotation=90,
        )

    output_path = output_path.expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        figure.tight_layout()
        figure.savefig(output_path, format="pdf", bbox_inches="tight")
    finally:
        plt.close(figure)


def plot_offline_random(input_dir: Path, output_path: Path) -> None:
    plot_results(load_archived_results(input_dir), output_path)


def plot_generated_offline_random(
    input_path: Path, output_path: Path
) -> pd.DataFrame:
    data = load_generated_results(input_path)
    plot_results(data, output_path)
    return data


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.input_dir is None:
        generated_input = args.generated_input or DEFAULT_GENERATED_INPUT
        data = plot_generated_offline_random(generated_input, args.out)
        incomplete = data[(data["eps"] <= data["d"]) & ~data["complete"]]
        if not incomplete.empty:
            details = ", ".join(
                f"{row.formula}(d={row.d}, eps={row.eps}: "
                f"{row.paired_samples}/{SAMPLES_PER_CELL})"
                for row in incomplete.itertuples()
            )
            print(f"omitted incomplete cells: {details}")
    else:
        plot_offline_random(args.input_dir, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
