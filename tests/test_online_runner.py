#!/usr/bin/env python3
"""Focused tests for the resumable online experiment runner."""

from __future__ import annotations

import contextlib
import csv
import importlib.util
import os
import io
import sys
import unittest
from collections import Counter
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import mock


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RUNNER_PATH = PROJECT_ROOT / "scripts/run_online.py"


def load_runner():
    specification = importlib.util.spec_from_file_location(
        "sttt_online_runner", RUNNER_PATH
    )
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load online runner from {RUNNER_PATH}")
    runner = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = runner
    specification.loader.exec_module(runner)
    return runner


runner = load_runner()


def make_job(
    job_id: str,
    *,
    mode: str = "incremental",
    formula: str = "ac",
    duration: int = 16,
    epsilon: int = 2,
    chunk_size: int = 4,
    sample: int = 0,
):
    stage = f"{mode}/{formula}"
    filename = runner.output_name(mode, formula)
    return runner.Job(
        job_id=job_id,
        stage=stage,
        label=job_id,
        target=f"sttt_online_{mode}_{formula}",
        arguments=("--repeat", "10"),
        output_path=Path(mode) / filename,
        temporary_path=Path(mode) / filename,
        mode=mode,
        formula=formula,
        duration=duration,
        epsilon=epsilon,
        chunk_size=chunk_size,
        sample=sample,
    )


def result_row(job, timing: float = 0.25, verdicts: str | None = None) -> str:
    if verdicts is None:
        verdicts = "2" * (job.duration // job.chunk_size)
    return (
        f"{job.duration} {job.epsilon} {job.chunk_size} {job.sample} "
        f"{timing} {verdicts}"
    )


class OnlinePlanTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.args = runner.parse_args([])
        cls.plan = runner.build_plan(cls.args)

    def test_default_plan_has_every_unique_configuration(self) -> None:
        counts = Counter(job.stage for job in self.plan)
        self.assertEqual(set(counts.values()), {3600})
        self.assertEqual(len(counts), 12)
        self.assertEqual(len(self.plan), 43200)
        self.assertEqual(len({job.job_id for job in self.plan}), len(self.plan))

    def test_default_plan_uses_ten_repetitions(self) -> None:
        for job in self.plan:
            with self.subTest(job=job.job_id):
                index = job.arguments.index("--repeat")
                self.assertEqual(job.arguments[index + 1], "10")

    def test_mode_and_formula_selection_is_deterministic(self) -> None:
        args = runner.parse_args(
            ["--mode", "naive", "--formula", "ad", "--formula", "ac"]
        )
        self.assertEqual(
            runner.selected_stages(args),
            ("naive/ac", "naive/ad"),
        )

    def test_command_passes_each_argument_once(self) -> None:
        job = self.plan[0]
        command = runner.command_for_job(job, Path("/tmp/build"), Path("/tmp/out"))
        for option in (
            "--duration",
            "--epsilon",
            "--chunk",
            "--sample",
            "--repeat",
            "--data-dir",
            "--output-dir",
            "--quiet",
        ):
            self.assertEqual(command.count(option), 1)
        output_index = command.index("--output-dir")
        self.assertEqual(command[output_index + 1], "/tmp/out/incremental")

    def test_config_tracks_build_context_and_no_build(self) -> None:
        args = runner.parse_args(
            ["--build-dir", "/tmp/sttt-online-build", "--no-build"]
        )
        context = {
            "build_directory": str(args.build_dir),
            "build_profile": "Release/O3",
            "cxx_compiler": "/toolchain/c++",
            "cxx_compiler_version": "test compiler",
            "cxxflags": "-march=native",
        }
        with mock.patch.object(
            runner, "build_context", return_value=context
        ) as context_for:
            config = runner.experiment_config(args)
        context_for.assert_called_once_with(args.build_dir)
        self.assertEqual(config["build"], context)
        self.assertTrue(config["no_build"])

    def test_build_selects_release_configuration(self) -> None:
        with mock.patch.object(runner.subprocess, "run") as run_process:
            runner.configure_and_build(
                ("incremental/ac",), Path("/tmp/sttt-online-build"), 2, False
            )
        self.assertEqual(run_process.call_count, 2)
        build_command = run_process.call_args_list[1].args[0]
        config_index = build_command.index("--config")
        self.assertEqual(build_command[config_index + 1], "Release")


class BuildContextTests(unittest.TestCase):
    def test_context_is_stable_when_cmake_cache_appears(self) -> None:
        with TemporaryDirectory() as temporary:
            build_directory = Path(temporary) / "build"
            compiler = Path(temporary) / "toolchain/c++"
            with mock.patch.dict(
                os.environ,
                {"CXX": str(compiler), "CXXFLAGS": "-march=native"},
            ):
                before = runner.build_context(build_directory)

            build_directory.mkdir()
            (build_directory / "CMakeCache.txt").write_text(
                f"CMAKE_CXX_COMPILER:FILEPATH={compiler}\n"
                "CMAKE_CXX_FLAGS:STRING=-march=native\n",
                encoding="utf-8",
            )
            with mock.patch.dict(
                os.environ,
                {"CXX": "/other/c++", "CXXFLAGS": "-O0"},
            ):
                after = runner.build_context(build_directory)

        self.assertEqual(before, after)
        self.assertEqual(before["build_profile"], "Release/O3")
        self.assertEqual(before["cxxflags"], "-march=native")


class ResultTests(unittest.TestCase):
    def test_validate_result_accepts_one_well_formed_row(self) -> None:
        job = make_job("incremental|ac|16|2|4|0")
        with TemporaryDirectory() as temporary:
            path = Path(temporary) / job.temporary_path
            path.parent.mkdir(parents=True)
            path.write_text(result_row(job) + "\n", encoding="utf-8")
            self.assertEqual(runner.validate_result(job, path), result_row(job))

    def test_validate_result_rejects_wrong_verdict_length(self) -> None:
        job = make_job("incremental|ac|16|2|4|0")
        with TemporaryDirectory() as temporary:
            path = Path(temporary) / job.temporary_path
            path.parent.mkdir(parents=True)
            path.write_text(result_row(job, verdicts="22") + "\n", encoding="utf-8")
            with self.assertRaisesRegex(RuntimeError, "verdict sequence"):
                runner.validate_result(job, path)

    def test_completed_timing_uses_per_repetition_average_seconds(self) -> None:
        job = make_job("incremental|ac|16|2|4|0")
        self.assertEqual(
            runner.format_attempt_timing(
                "complete", result_row(job), 99.0, (4,), "s"
            ),
            "avg 0.25s",
        )

    def test_summary_uses_mean_of_sample_speedups(self) -> None:
        plan = []
        records = {}
        for sample in runner.SAMPLES:
            incremental = make_job(
                f"incremental|ac|64|2|4|{sample}",
                duration=64,
                sample=sample,
            )
            naive = make_job(
                f"naive|ac|64|2|4|{sample}",
                mode="naive",
                duration=64,
                sample=sample,
            )
            plan.extend((incremental, naive))
            records[incremental.job_id] = {
                "status": "complete",
                "result_row": result_row(incremental, timing=2.0),
            }
            records[naive.job_id] = {
                "status": "complete",
                "result_row": result_row(naive, timing=4.0),
            }

        with TemporaryDirectory() as temporary:
            output = Path(temporary)
            runner.materialize_summary(plan, records, output)
            with (output / runner.SUMMARY_FILENAME).open(
                encoding="utf-8"
            ) as source:
                rows = list(csv.DictReader(source))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["formula"], "Psi 1")
        self.assertEqual(float(rows[0]["speedup"]), 2.0)
        self.assertEqual(float(rows[0]["incr wall time"]), 2.0)
        self.assertEqual(float(rows[0]["time per portion"]), 0.125)


class CheckpointTests(unittest.TestCase):
    def test_main_reports_sqlite_errors(self) -> None:
        with mock.patch.object(
            runner,
            "run",
            side_effect=runner.sqlite3.OperationalError("broken checkpoint"),
        ), contextlib.redirect_stderr(io.StringIO()) as errors:
            self.assertEqual(runner.main([]), 1)
        self.assertIn("error: broken checkpoint", errors.getvalue())

    def test_interrupted_run_resumes_at_next_configuration(self) -> None:
        jobs = [
            make_job("incremental|ac|16|2|4|0", sample=0),
            make_job("incremental|ac|16|2|4|1", sample=1),
        ]
        rows = tuple(result_row(job) for job in jobs)
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            args = runner.parse_args(
                [
                    "--mode",
                    "incremental",
                    "--formula",
                    "ac",
                    "--output-dir",
                    str(root / "results"),
                    "--build-dir",
                    str(root / "build"),
                    "--no-build",
                ]
            )
            fingerprint = mock.patch.object(
                runner,
                "experiment_fingerprint",
                return_value=("fingerprint", "source", "config"),
            )
            with mock.patch.object(runner, "build_plan", return_value=jobs), \
                    mock.patch.object(runner, "configure_and_build"), \
                    fingerprint, \
                    mock.patch.object(
                        runner,
                        "execute_job",
                        side_effect=[
                            runner.Attempt("complete", rows[0], 0.1, ""),
                            KeyboardInterrupt(),
                        ],
                    ) as execute:
                with contextlib.redirect_stdout(io.StringIO()):
                    self.assertEqual(runner.run(args), 130)
                self.assertEqual(execute.call_count, 2)

            with mock.patch.object(runner, "build_plan", return_value=jobs), \
                    mock.patch.object(runner, "configure_and_build"), \
                    mock.patch.object(
                        runner,
                        "experiment_fingerprint",
                        return_value=("fingerprint", "source", "config"),
                    ), \
                    mock.patch.object(
                        runner,
                        "execute_job",
                        return_value=runner.Attempt("complete", rows[1], 0.2, ""),
                    ) as execute:
                with contextlib.redirect_stdout(io.StringIO()):
                    self.assertEqual(runner.run(args), 0)
                execute.assert_called_once()
                self.assertEqual(execute.call_args.args[0], jobs[1])

            result = args.output_dir / jobs[0].output_path
            self.assertEqual(
                result.read_text(encoding="utf-8"),
                "".join(f"{row}\n" for row in rows),
            )


if __name__ == "__main__":
    unittest.main()
