#!/usr/bin/env python3
"""Regression tests for the offline random-trace heatmap aggregation."""

from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import mock


PROJECT_ROOT = Path(__file__).resolve().parents[1]
PLOTTER_PATH = PROJECT_ROOT / "visualization/plot_offline_random.py"


def load_plotter():
    specification = importlib.util.spec_from_file_location(
        "sttt_plot_offline_random", PLOTTER_PATH
    )
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load heatmap plotter from {PLOTTER_PATH}")
    plotter = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(plotter)
    return plotter


plotter = load_plotter()


class OfflineRandomAggregationTests(unittest.TestCase):
    def test_verified_timeouts_affect_both_speedups_but_not_false_positives(
        self,
    ) -> None:
        with TemporaryDirectory() as temporary_directory:
            result_root = Path(temporary_directory)
            approximate_path = result_root / "approximate/random/phi6.txt"
            exact_path = result_root / "exact/random/phi6.txt"
            approximate_path.parent.mkdir(parents=True)
            exact_path.parent.mkdir(parents=True)

            timeout_verdicts = {
                (1, 12): 0,
                (1, 35): 2,
                (1, 48): 0,
                (1, 79): 0,
                (1, 97): 2,
                (2, 48): 2,
            }
            approximate_rows = []
            exact_rows = []
            for epsilon in (1, 2):
                for sample in range(100):
                    approximate_verdict = timeout_verdicts.get(
                        (epsilon, sample), 0
                    )
                    if epsilon == 1 and sample == 0:
                        approximate_verdict = 2
                    approximate_rows.append(
                        f"8 {epsilon} 0 {sample} 0 1 {approximate_verdict}"
                    )
                    if (epsilon, sample) in timeout_verdicts:
                        continue
                    exact_rows.append(
                        f"8 {epsilon} - {sample} - 2 1 0 1 0"
                    )

            approximate_path.write_text(
                "\n".join(approximate_rows) + "\n", encoding="utf-8"
            )
            exact_path.write_text(
                "\n".join(exact_rows) + "\n", encoding="utf-8"
            )
            timeout_rows = [
                "stage\tinstance\ttimeout_seconds\treason",
                *(
                    "exact-random\t"
                    f"exact-random|phi6|8|{epsilon}|{sample}|grid-2d-v1"
                    "\t600\tZ3 timed out: timeout"
                    for epsilon, sample in timeout_verdicts
                ),
            ]
            (result_root / "timeouts.tsv").write_text(
                "\n".join(timeout_rows) + "\n", encoding="utf-8"
            )

            with (
                mock.patch.object(plotter, "FORMULAS", ("phi6",)),
                mock.patch.object(plotter, "DURATIONS", (8,)),
                mock.patch.object(plotter, "EPSILON_VALUES", (1, 2)),
            ):
                data = plotter.load_generated_results(result_root)

        epsilon_1 = data[data["eps"] == 1].iloc[0]
        self.assertTrue(epsilon_1["complete"])
        self.assertEqual(epsilon_1["paired_samples"], 100)
        self.assertEqual(epsilon_1["imputed_timeouts"], 5)
        self.assertAlmostEqual(epsilon_1["ADM"], 1.0)
        self.assertAlmostEqual(epsilon_1["EDM"], 31.9)
        self.assertAlmostEqual(epsilon_1["FP"], 1.0)
        self.assertAlmostEqual(epsilon_1["speedup"], 31.9)
        self.assertAlmostEqual(epsilon_1["combined"], 3190.0 / 1302.0)

        epsilon_2 = data[data["eps"] == 2].iloc[0]
        self.assertTrue(epsilon_2["complete"])
        self.assertEqual(epsilon_2["paired_samples"], 100)
        self.assertEqual(epsilon_2["imputed_timeouts"], 1)
        self.assertAlmostEqual(epsilon_2["ADM"], 1.0)
        self.assertAlmostEqual(epsilon_2["EDM"], 7.98)
        self.assertAlmostEqual(epsilon_2["FP"], 0.0)
        self.assertAlmostEqual(epsilon_2["speedup"], 7.98)
        self.assertAlmostEqual(epsilon_2["combined"], 798.0 / 700.0)


if __name__ == "__main__":
    unittest.main()
