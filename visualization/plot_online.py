#!/usr/bin/env python3
"""Generate heatmaps from the online benchmark results."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm, Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable


BENCHMARK_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = BENCHMARK_ROOT / "results/reference/online/summary.csv"
DEFAULT_OUTPUT = BENCHMARK_ROOT / "results/figures/speedup_online.pdf"


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.rename(columns={column: column.strip().lower() for column in df.columns})
    aliases = {
        "incr_wall_time": "incr wall time",
        "time_per_portion": "time per portion",
    }
    for source, destination in aliases.items():
        if source in df.columns and destination not in df.columns:
            df = df.rename(columns={source: destination})

    required = [
        "d",
        "eps",
        "delta",
        "speedup",
        "incr wall time",
        "time per portion",
    ]
    for column in required:
        if column not in df.columns:
            raise ValueError(f"Missing column: '{column}'")

    if "formula" not in df.columns:
        df["formula"] = "Psi"

    for column in required:
        df[column] = pd.to_numeric(df[column], errors="coerce")

    return df


def build_grid(df: pd.DataFrame, value_column: str):
    durations = sorted(df["d"].dropna().unique().tolist())
    pairs_df = (
        df[["eps", "delta"]]
        .dropna()
        .drop_duplicates()
        .sort_values(["eps", "delta"])
    )
    pairs = pairs_df[["eps", "delta"]].values.tolist()
    values = np.full((len(pairs), len(durations)), np.nan)

    for row, (epsilon, delta) in enumerate(pairs):
        for column, duration in enumerate(durations):
            selected = df[
                (df["eps"] == epsilon)
                & (df["delta"] == delta)
                & (df["d"] == duration)
            ]
            if len(selected) == 1:
                values[row, column] = float(selected[value_column].values[0])

    y_labels = [f"{int(epsilon)}/{int(delta)}" for epsilon, delta in pairs]
    x_labels = [str(int(duration)) for duration in durations]
    return values, y_labels, x_labels


def format_cell(metric: str, value: float) -> str:
    if not np.isfinite(value):
        return ""
    if metric == "speedup":
        return f"{value:.2f}"

    milliseconds = value * 1e3
    return f"{milliseconds:.2f}" if milliseconds < 10 else f"{milliseconds:.1f}"


def text_color(cmap, norm, value: float) -> str:
    value = max(value, np.finfo(float).eps)
    red, green, blue, _ = cmap(norm(value))
    luminance = 0.2126 * red + 0.7152 * green + 0.0722 * blue
    return "black" if luminance > 0.6 else "white"


def draw_heatmap(
    axis,
    values,
    x_labels,
    y_labels,
    cmap_name,
    log_scale=True,
):
    minimum = np.nanmin(values)
    maximum = np.nanmax(values)
    if log_scale:
        minimum = max(minimum, np.finfo(float).eps)
    norm = (
        LogNorm(vmin=minimum, vmax=maximum)
        if log_scale
        else Normalize(vmin=minimum, vmax=maximum)
    )
    cmap = plt.get_cmap(cmap_name)
    image = axis.imshow(values, origin="lower", aspect="auto", cmap=cmap, norm=norm)

    for spine in axis.spines.values():
        spine.set_visible(False)

    rows, columns = values.shape
    axis.set_xticks(np.arange(-0.5, columns, 1), minor=True)
    axis.set_yticks(np.arange(-0.5, rows, 1), minor=True)
    axis.grid(which="minor", color="white", linestyle="-", linewidth=0.6)
    axis.tick_params(which="minor", bottom=False, left=False)

    axis.set_xticks(range(columns))
    axis.set_xticklabels(x_labels, fontsize=9)
    axis.set_yticks(range(rows))
    axis.set_yticklabels(y_labels, fontsize=9)

    divider = make_axes_locatable(axis)
    colorbar_axis = divider.append_axes("right", size="3.5%", pad=0.05)
    colorbar = plt.colorbar(image, cax=colorbar_axis)
    colorbar.ax.tick_params(labelsize=8)

    return cmap, norm


def plot_grid(
    df: pd.DataFrame,
    output_path: Path,
    cmap_name="GnBu",
    log_scale=True,
    figure_width=10.5,
    row_height=2.0,
    keep_d16=False,
):
    excluded_durations = [32] if keep_d16 else [16, 32]
    df = df[~df["d"].isin(excluded_durations)].copy()

    formulas = list(df["formula"].dropna().unique())
    metrics = [
        ("speedup", "Speed-up (times)"),
        ("incr wall time", "Total time (ms)"),
        ("time per portion", "Per-portion time (ms)"),
    ]
    row_count = len(formulas)
    column_count = len(metrics)

    figure = plt.figure(figsize=(figure_width, row_height * row_count))
    plt.subplots_adjust(
        left=0.30,
        right=0.985,
        top=0.95,
        bottom=0.08,
        wspace=0.16,
        hspace=0.35,
    )

    axes = np.empty((row_count, column_count), dtype=object)
    for row, formula in enumerate(formulas):
        subset = df[df["formula"] == formula]
        speedup, y_labels, x_labels = build_grid(subset, "speedup")
        incremental_time, _, _ = build_grid(subset, "incr wall time")
        portion_time, _, _ = build_grid(subset, "time per portion")
        metric_values = [speedup, incremental_time, portion_time]

        for column, (metric, _) in enumerate(metrics):
            axis = plt.subplot(
                row_count,
                column_count,
                row * column_count + column + 1,
            )
            axes[row, column] = axis
            cmap, norm = draw_heatmap(
                axis,
                metric_values[column],
                x_labels=x_labels,
                y_labels=y_labels if column == 0 else [""] * len(y_labels),
                cmap_name=cmap_name,
                log_scale=log_scale,
            )

            values = metric_values[column]
            for value_row in range(values.shape[0]):
                for value_column in range(values.shape[1]):
                    value = values[value_row, value_column]
                    if not np.isfinite(value):
                        continue
                    axis.text(
                        value_column,
                        value_row,
                        format_cell(metric, value),
                        ha="center",
                        va="center",
                        fontsize=8,
                        color=text_color(cmap, norm, value),
                    )

            axis.set_xlabel("duration", fontsize=10, labelpad=4)
            if column == 0:
                axis.set_ylabel("epsilon/Delta", fontsize=10)
                axis.yaxis.set_label_coords(-0.18, 0.5)
                axis.text(
                    -0.36,
                    0.5,
                    str(formula),
                    rotation=90,
                    transform=axis.transAxes,
                    ha="right",
                    va="center",
                    fontsize=11,
                )
            else:
                axis.set_yticklabels([])

    for column, (_, title) in enumerate(metrics):
        first_axis = axes[0, column]
        left, _, width, _ = first_axis.get_position().bounds
        figure.text(
            left + width / 2.0,
            0.975,
            title,
            ha="center",
            va="top",
            fontsize=12,
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    if output_path.suffix.lower() != ".pdf":
        figure.savefig(output_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(figure)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot heatmaps from the online benchmark results."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help=f"input CSV or TSV file (default: {DEFAULT_INPUT})",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"output figure path (default: {DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "--delimiter",
        default=None,
        help="override delimiter (default: .tsv -> tab, otherwise comma)",
    )
    parser.add_argument("--cmap", default="GnBu")
    parser.add_argument("--no-log", action="store_true")
    parser.add_argument("--fig-width", type=float, default=10.5)
    parser.add_argument("--row-height", type=float, default=2.0)
    parser.add_argument(
        "--keep-d16",
        action="store_true",
        help="keep d=16; d=32 is always omitted",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    input_path = args.input.expanduser().resolve()
    output_path = args.out.expanduser().resolve()
    delimiter = (
        args.delimiter
        if args.delimiter is not None
        else "\t" if input_path.suffix.lower() == ".tsv" else ","
    )
    data = normalize_columns(pd.read_csv(input_path, sep=delimiter))
    plot_grid(
        df=data,
        output_path=output_path,
        cmap_name=args.cmap,
        log_scale=not args.no_log,
        figure_width=args.fig_width,
        row_height=args.row_height,
        keep_d16=args.keep_d16,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
