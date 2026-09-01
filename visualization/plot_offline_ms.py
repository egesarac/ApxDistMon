#!/usr/bin/env python3
"""Generate the offline mutual-separation case-study figure."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


BENCHMARK_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_GENERATED_INPUT = BENCHMARK_ROOT / "results/reference/offline"
DEFAULT_OUTPUT = BENCHMARK_ROOT / "results/figures/speedup_offline_ms.pdf"
EPSILON_VALUES = (0.05, 0.1, 0.15, 0.2, 0.25)
AGENT_COUNTS = (2, 3, 4)
METHODS = (
    ("EDM", "SMT", "o"),
    ("ADM", "ORIG", "s"),
    ("ADM-F", "FINE", "D"),
    ("ADM-FR", "F-RELATIVE", "^"),
)
GENERATED_FILES = (
    ("ADM", "mutual_separation.txt"),
    ("ADM-F", "mutual_separation_adm_f.txt"),
    ("ADM-FR", "mutual_separation_adm_fr.txt"),
)
SAMPLING_PERIOD_SECONDS = 0.05
EXACT_TIMEOUT_SECONDS = 360.0


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group()
    source.add_argument(
        "--input",
        type=Path,
        help="archived result workbook",
    )
    source.add_argument(
        "--generated-input",
        type=Path,
        help=f"result root (default: {DEFAULT_GENERATED_INPUT})",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"output figure path (default: {DEFAULT_OUTPUT})",
    )
    return parser.parse_args(argv)


def _row_index(table: pd.DataFrame, label: str) -> int:
    labels = table.iloc[:, 0].fillna("").astype(str).str.strip()
    matches = table.index[labels == label].tolist()
    if len(matches) != 1:
        raise ValueError(f"expected one {label!r} row, found {len(matches)}")
    return int(matches[0])


def load_results(input_path: Path) -> dict[int, dict[str, list[float]]]:
    input_path = input_path.expanduser().resolve()
    if not input_path.is_file():
        raise FileNotFoundError(f"missing mutual-separation results: {input_path}")

    table = pd.read_excel(input_path, header=None)
    agents_row = _row_index(table, "n =")
    epsilon_row = _row_index(table, "eps (in seconds) =")
    method_rows = {
        method: _row_index(table, source_label)
        for method, source_label, _ in METHODS
    }

    results: dict[int, dict[str, list[float]]] = {}
    block_width = len(EPSILON_VALUES)
    for block, agents in enumerate(AGENT_COUNTS):
        start = 1 + block * block_width
        observed_agents = int(table.iloc[agents_row, start])
        if observed_agents != agents:
            raise ValueError(
                f"agent block {block} is n={observed_agents}; expected n={agents}"
            )

        observed_epsilons = tuple(
            float(value)
            for value in table.iloc[
                epsilon_row, start : start + block_width
            ]
        )
        if observed_epsilons != EPSILON_VALUES:
            raise ValueError(
                f"n={agents} epsilon values are {observed_epsilons}; "
                f"expected {EPSILON_VALUES}"
            )

        results[agents] = {}
        for method, _, _ in METHODS:
            values: list[float] = []
            for raw_value in table.iloc[
                method_rows[method], start : start + block_width
            ]:
                if str(raw_value).strip().lower() == "oot":
                    values.append(float("nan"))
                    continue
                value = float(raw_value)
                if not math.isfinite(value) or value <= 0:
                    raise ValueError(
                        f"{method}, n={agents} contains invalid time {raw_value!r}"
                    )
                values.append(value)
            results[agents][method] = values
    return results


def _result_rows(path: Path) -> list[list[str]]:
    if not path.is_file():
        raise FileNotFoundError(f"missing generated result: {path}")
    return [
        line.split()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]


def load_generated_results(
    input_path: Path,
) -> tuple[
    tuple[float, ...],
    tuple[int, ...],
    dict[int, dict[str, list[float]]],
]:
    input_path = input_path.expanduser().resolve()
    official_approximate = input_path / "approximate/case_studies"
    if official_approximate.is_dir():
        run_paths = [(
            official_approximate,
            input_path / "exact/case_studies/mutual_separation.txt",
        )]
    else:
        run_paths = [
            (
                approximate_directory,
                approximate_directory.parent /
                "exact/case_studies/mutual_separation.txt",
            )
            for approximate_directory in sorted(
                input_path.glob("**/approximate")
            )
        ]
    if not run_paths:
        raise FileNotFoundError(
            f"no generated mutual-separation runs in {input_path}"
        )

    observations: dict[tuple[int, int, str], float] = {}
    for approximate_directory, exact_path in run_paths:
        for method, filename in GENERATED_FILES:
            for fields in _result_rows(approximate_directory / filename):
                agents = int(fields[1])
                multiplier = int(fields[2])
                observations[(agents, multiplier, method)] = float(fields[4])

        for fields in _result_rows(exact_path):
            agents = int(fields[1])
            multiplier = int(fields[2])
            observations[(agents, multiplier, "EDM")] = max(
                float(fields[3]),
                float(fields[5]),
            )

    agent_counts = tuple(sorted({agents for agents, _, _ in observations}))
    multipliers = tuple(sorted({multiplier for _, multiplier, _ in observations}))
    epsilon_values = tuple(
        round(multiplier * SAMPLING_PERIOD_SECONDS, 10)
        for multiplier in multipliers
    )
    results = {
        agents: {
            method: [
                observations.get((agents, multiplier, method), float("nan"))
                for multiplier in multipliers
            ]
            for method, _, _ in METHODS
        }
        for agents in agent_counts
    }
    return epsilon_values, agent_counts, results


def plot_results(
    results: dict[int, dict[str, list[float]]],
    epsilon_values: tuple[float, ...],
    agent_counts: tuple[int, ...],
    output_path: Path,
) -> None:
    palette = sns.color_palette("colorblind", n_colors=len(METHODS))
    figure, axes_grid = plt.subplots(
        1,
        len(agent_counts),
        figsize=(4 * len(agent_counts), 4),
        squeeze=False,
    )
    axes = axes_grid[0]

    for index, agents in enumerate(agent_counts):
        axis = axes[index]
        for (method, _, marker), color in zip(METHODS, palette):
            axis.plot(
                epsilon_values,
                np.asarray(results[agents][method], dtype=float),
                color=color,
                label=method,
                marker=marker,
                linestyle="-",
            )

        axis.set_xlabel("ε (seconds)", fontsize=12)
        if index == 0:
            axis.set_ylabel("Time (seconds)", fontsize=12)
        axis.set_title(f"n={agents}", fontsize=14)
        axis.set_yscale("log")
        axis.set_ylim(top=1e3)
        axis.set_xlim(
            min(epsilon_values) * 0.9,
            max(epsilon_values) * 1.1,
        )
        axis.set_xticks(epsilon_values)
        axis.set_xticklabels(epsilon_values)
        axis.tick_params(axis="x", labelsize=9)
        axis.grid(True, which="both", linestyle="--", linewidth=0.5)

        exact_times = np.asarray(results[agents]["EDM"], dtype=float)
        timed_out = ~np.isfinite(exact_times)
        if np.any(timed_out):
            axis.scatter(
                np.asarray(epsilon_values)[timed_out],
                np.full(np.count_nonzero(timed_out), EXACT_TIMEOUT_SECONDS),
                color=palette[0],
                marker="X",
                s=50,
                zorder=4,
            )

    handles, labels = axes[0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        loc="upper center",
        ncol=len(METHODS),
        fontsize=8,
        frameon=False,
    )

    output_path = output_path.expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        figure.tight_layout(rect=(0.035, 0, 1, 0.88))
        figure.savefig(output_path, bbox_inches="tight")
    finally:
        plt.close(figure)


def plot_offline_ms(input_path: Path, output_path: Path) -> None:
    plot_results(
        load_results(input_path),
        EPSILON_VALUES,
        AGENT_COUNTS,
        output_path,
    )


def plot_generated_offline_ms(input_path: Path, output_path: Path) -> None:
    epsilon_values, agent_counts, results = load_generated_results(input_path)
    plot_results(results, epsilon_values, agent_counts, output_path)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.input is not None:
        plot_offline_ms(args.input, args.out)
    else:
        plot_generated_offline_ms(
            args.generated_input or DEFAULT_GENERATED_INPUT,
            args.out,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
