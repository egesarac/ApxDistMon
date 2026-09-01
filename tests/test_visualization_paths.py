#!/usr/bin/env python3
"""Regression tests for layout-independent visualization defaults."""

from __future__ import annotations

import importlib.util
import shutil
import unittest
from itertools import count
from pathlib import Path
from tempfile import TemporaryDirectory


BENCHMARK_ROOT = Path(__file__).resolve().parents[1]
VISUALIZATION_ROOT = BENCHMARK_ROOT / "visualization"
MODULE_NAMES = (
    "plot_online",
    "plot_offline_random",
    "plot_offline_ms",
    "plot_offline_wt",
)
EXPECTED_DEFAULTS = {
    "plot_online": {
        "DEFAULT_INPUT": "results/reference/online/summary.csv",
        "DEFAULT_OUTPUT": "results/figures/speedup_online.pdf",
    },
    "plot_offline_random": {
        "DEFAULT_GENERATED_INPUT": "results/reference/offline",
        "DEFAULT_OUTPUT": "results/figures/speedup_offline_random.pdf",
    },
    "plot_offline_ms": {
        "DEFAULT_GENERATED_INPUT": "results/reference/offline",
        "DEFAULT_OUTPUT": "results/figures/speedup_offline_ms.pdf",
    },
    "plot_offline_wt": {
        "DEFAULT_GENERATED_INPUT": "results/reference/offline",
        "DEFAULT_OUTPUT": "results/figures/speedup_offline_wt.pdf",
    },
}
EXPLICIT_INPUT_OPTIONS = {
    "plot_online": ("--input", "input"),
    "plot_offline_random": ("--input-dir", "input_dir"),
    "plot_offline_ms": ("--input", "input"),
    "plot_offline_wt": ("--input", "input"),
}
_IMPORT_SEQUENCE = count()


def load_module(module_name: str, visualization_root: Path):
    module_path = visualization_root / f"{module_name}.py"
    specification = importlib.util.spec_from_file_location(
        f"visualization_path_test_{module_name}_{next(_IMPORT_SEQUENCE)}",
        module_path,
    )
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load visualization module from {module_path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


class VisualizationPathTests(unittest.TestCase):
    def assert_module_paths(
        self,
        module_name: str,
        visualization_root: Path,
        expected_root: Path,
    ) -> None:
        module = load_module(module_name, visualization_root)
        self.assertEqual(module.BENCHMARK_ROOT, expected_root.resolve())
        self.assertFalse(hasattr(module, "REPOSITORY_ROOT"))

        for constant, relative_path in EXPECTED_DEFAULTS[module_name].items():
            self.assertEqual(
                getattr(module, constant),
                expected_root.resolve() / relative_path,
            )

        default_arguments = module.parse_args([])
        self.assertEqual(default_arguments.out, module.DEFAULT_OUTPUT)
        if module_name == "plot_online":
            self.assertEqual(default_arguments.input, module.DEFAULT_INPUT)
        else:
            self.assertIsNone(default_arguments.generated_input)

        option, destination = EXPLICIT_INPUT_OPTIONS[module_name]
        legacy_path = Path("../legacy-results/input")
        arguments = module.parse_args([option, str(legacy_path)])
        self.assertEqual(getattr(arguments, destination), legacy_path)

        if module_name == "plot_offline_random":
            self.assertFalse(hasattr(module, "DEFAULT_INPUT_DIR"))
        elif module_name in {"plot_offline_ms", "plot_offline_wt"}:
            self.assertFalse(hasattr(module, "DEFAULT_INPUT"))

    def test_defaults_are_relative_to_current_benchmark_root(self) -> None:
        for module_name in MODULE_NAMES:
            with self.subTest(module=module_name):
                self.assert_module_paths(
                    module_name,
                    VISUALIZATION_ROOT,
                    BENCHMARK_ROOT,
                )

    def test_defaults_are_relative_to_standalone_repository_root(self) -> None:
        with TemporaryDirectory() as temporary_directory:
            standalone_root = Path(temporary_directory)
            standalone_visualization = standalone_root / "visualization"
            standalone_visualization.mkdir()
            for module_name in MODULE_NAMES:
                shutil.copy2(
                    VISUALIZATION_ROOT / f"{module_name}.py",
                    standalone_visualization,
                )

            for module_name in MODULE_NAMES:
                with self.subTest(module=module_name):
                    self.assert_module_paths(
                        module_name,
                        standalone_visualization,
                        standalone_root,
                    )


if __name__ == "__main__":
    unittest.main()
