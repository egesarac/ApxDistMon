#!/usr/bin/env python3
"""Focused tests for the resumable offline experiment runner."""

from __future__ import annotations

import contextlib
import importlib.util
import io
import json
import sys
import unittest
from collections import Counter
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import mock


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RUNNER_PATH = PROJECT_ROOT / "scripts/run_offline.py"


def load_runner():
    specification = importlib.util.spec_from_file_location(
        "sttt_offline_runner", RUNNER_PATH
    )
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load offline runner from {RUNNER_PATH}")
    runner = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = runner
    specification.loader.exec_module(runner)
    return runner


runner = load_runner()


def make_job(
    job_id: str,
    *,
    stage: str = "approximate-random",
    output_path: Path = Path("approximate/random/phi1.txt"),
    timeout_seconds: int | None = None,
):
    return runner.Job(
        job_id=job_id,
        stage=stage,
        label=job_id,
        entrypoint="unused",
        arguments=("--repeat", "30"),
        output_path=output_path,
        temporary_path=Path("phi1.txt"),
        timeout_seconds=timeout_seconds,
    )


class OfflinePlanTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.args = runner.parse_args([])
        cls.plan = runner.build_plan(cls.args)

    def test_default_plan_has_one_unique_job_per_configuration(self) -> None:
        expected = {
            "approximate-random": 9000,
            "approximate-mutual-separation": 45,
            "approximate-water-tanks": 60,
            "exact-random": 7400,
            "exact-mutual-separation": 15,
            "exact-water-tanks": 12,
        }
        self.assertEqual(Counter(job.stage for job in self.plan), expected)
        self.assertEqual(len(self.plan), 16532)
        self.assertEqual(len({job.job_id for job in self.plan}), len(self.plan))

    def test_random_plan_uses_reported_duration_and_epsilon_matrix(self) -> None:
        groups = (
            (
                "approximate",
                [
                    job
                    for job in self.plan
                    if job.stage == "approximate-random"
                ],
                {4, 8, 16, 32},
            ),
            (
                "exact-untimed",
                [
                    job
                    for job in self.plan
                    if job.stage == "exact-random"
                    and job.job_id.split("|")[1]
                    not in runner.EXACT_TIMED_RANDOM_FORMULAS
                ],
                {4, 8, 16, 32},
            ),
            (
                "exact-timed",
                [
                    job
                    for job in self.plan
                    if job.stage == "exact-random"
                    and job.job_id.split("|")[1]
                    in runner.EXACT_TIMED_RANDOM_FORMULAS
                ],
                {4, 8},
            ),
        )
        for name, jobs, expected_durations in groups:
            durations = {
                int(job.job_id.split("|")[2])
                for job in jobs
            }
            epsilons = {
                int(job.job_id.split("|")[3])
                for job in jobs
            }
            with self.subTest(group=name):
                self.assertEqual(durations, expected_durations)
                self.assertEqual(epsilons, {1, 2, 4, 8})

    def test_case_studies_run_after_exact_random(self) -> None:
        ordered = runner.jobs_for_stages(self.plan, runner.STAGES)
        observed = tuple(dict.fromkeys(job.stage for job in ordered))
        self.assertEqual(observed, runner.STAGES)
        self.assertEqual(observed[:2], ("approximate-random", "exact-random"))

    def test_exact_case_study_job_ids_version_the_backend(self) -> None:
        jobs = [
            job
            for job in self.plan
            if job.stage in {
                "exact-mutual-separation",
                "exact-water-tanks",
            }
        ]
        self.assertTrue(
            all(
                job.job_id.endswith("|rv-full-flow-v1")
                for job in jobs
            )
        )

    def test_exact_timed_random_job_ids_version_the_backend(self) -> None:
        jobs = [
            job
            for job in self.plan
            if job.stage == "exact-random"
            and job.job_id.split("|")[1]
            in runner.EXACT_TIMED_RANDOM_FORMULAS
        ]
        self.assertTrue(
            all(job.job_id.endswith("|grid-2d-v1") for job in jobs)
        )

    def test_formula_and_duration_filters_select_random_jobs(self) -> None:
        args = runner.parse_args(
            [
                "--only",
                "exact-random",
                "--formula",
                "phi5",
                "--formula",
                "phi6",
                "--duration",
                "4",
                "--duration",
                "8",
            ]
        )
        selected = runner.select_jobs(
            self.plan, runner.selected_stages(args), args
        )
        self.assertEqual(len(selected), 1400)
        self.assertEqual(
            {job.job_id.split("|")[1] for job in selected},
            {"phi5", "phi6"},
        )
        self.assertEqual(
            {int(job.job_id.split("|")[2]) for job in selected},
            {4, 8},
        )

    def test_timeout_retry_is_included_outside_random_filters(self) -> None:
        args = runner.parse_args(
            [
                "--only",
                "exact-random",
                "--formula",
                "phi5",
                "--duration",
                "4",
                "--rerun-timeouts",
            ]
        )
        timeout_job_id = "exact-random|phi2|32|8|37"
        selected = runner.select_jobs(
            self.plan,
            runner.selected_stages(args),
            args,
            {timeout_job_id},
        )
        self.assertIn(timeout_job_id, {job.job_id for job in selected})
        self.assertEqual(len(selected), 301)

    def test_default_repetitions_and_approximate_timeout_match_paper(self) -> None:
        expected_repetitions = {
            "approximate-random": "100",
            "approximate-mutual-separation": "10",
            "approximate-water-tanks": "10",
            "exact-random": "3",
            "exact-mutual-separation": "10",
            "exact-water-tanks": "10",
        }
        for job in self.plan:
            with self.subTest(job=job.job_id):
                repeat_index = job.arguments.index("--repeat")
                self.assertEqual(
                    job.arguments[repeat_index + 1],
                    expected_repetitions[job.stage],
                )
                if job.stage.startswith("approximate-"):
                    self.assertNotIn("--timeout", job.arguments)
                    self.assertIsNone(job.timeout_seconds)

    def test_case_study_option_covers_approximate_and_exact_jobs(self) -> None:
        args = runner.parse_args(
            [
                "--approximate-repetitions",
                "11",
                "--exact-repetitions",
                "13",
                "--case-study-repetitions",
                "17",
            ]
        )
        plan = runner.build_plan(args)
        for job in plan:
            expected = (
                "11"
                if job.stage == "approximate-random"
                else "13" if job.stage == "exact-random" else "17"
            )
            with self.subTest(job=job.job_id):
                repeat_index = job.arguments.index("--repeat")
                self.assertEqual(job.arguments[repeat_index + 1], expected)
        self.assertEqual(
            runner.experiment_config(args)["case_study_repetitions"],
            17,
        )

    def test_exact_jobs_pass_smt_timeout(self) -> None:
        expected_timeouts = {
            "exact-random": 600,
            "exact-mutual-separation": 360,
            "exact-water-tanks": 120,
        }
        for job in self.plan:
            if not job.stage.startswith("exact-"):
                continue
            with self.subTest(job=job.job_id):
                timeout_index = job.arguments.index("--timeout")
                expected = expected_timeouts[job.stage]
                self.assertEqual(job.arguments[timeout_index + 1], str(expected))
                self.assertEqual(job.timeout_seconds, expected)

    def test_config_tracks_build_context_and_no_build(self) -> None:
        args = runner.parse_args(
            ["--build-dir", "/tmp/sttt-offline-build", "--no-build"]
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


class CheckpointTests(unittest.TestCase):
    @staticmethod
    def timeout_config(random_timeout: int, exact_repetitions: int = 3) -> str:
        return json.dumps(
            {
                "exact_repetitions": exact_repetitions,
                "random_timeout": random_timeout,
                "mutual_separation_timeout": 360,
                "water_tanks_timeout": 120,
            },
            sort_keys=True,
            separators=(",", ":"),
        )

    def test_interrupted_run_resumes_at_the_next_configuration(self) -> None:
        jobs = [
            make_job("approximate-random|phi1|4|1|0"),
            make_job("approximate-random|phi1|4|1|1"),
        ]
        rows = (
            "4 1 0 0 7 0.1 1",
            "4 1 0 1 8 0.2 2",
        )
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            args = runner.parse_args(
                [
                    "--only",
                    "approximate-random",
                    "--output-dir",
                    str(root / "results"),
                    "--build-dir",
                    str(root / "build"),
                    "--no-build",
                ]
            )
            common_patches = (
                mock.patch.object(runner, "build_plan", return_value=jobs),
                mock.patch.object(runner, "configure_and_build"),
                mock.patch.object(
                    runner,
                    "experiment_fingerprint",
                    return_value=("fingerprint", "source", "config"),
                ),
            )

            with common_patches[0], common_patches[1], common_patches[2]:
                with mock.patch.object(
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

            state_path = args.output_dir / runner.STATE_FILENAME
            connection = runner.connect_checkpoint(state_path)
            try:
                records = runner.load_records(connection)
            finally:
                connection.close()
            self.assertEqual(set(records), {jobs[0].job_id})
            result_path = args.output_dir / jobs[0].output_path
            self.assertEqual(result_path.read_text(encoding="utf-8"), rows[0] + "\n")

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
                        return_value=runner.Attempt(
                            "complete", rows[1], 0.2, ""
                        ),
                    ) as execute:
                with contextlib.redirect_stdout(io.StringIO()):
                    self.assertEqual(runner.run(args), 0)
                execute.assert_called_once()
                self.assertEqual(execute.call_args.args[0], jobs[1])

            connection = runner.connect_checkpoint(state_path)
            try:
                records = runner.load_records(connection)
            finally:
                connection.close()
            self.assertEqual(set(records), {job.job_id for job in jobs})
            self.assertEqual(
                result_path.read_text(encoding="utf-8"),
                "".join(f"{row}\n" for row in rows),
            )

    def test_timeout_is_checkpointed_reported_and_omitted_from_results(self) -> None:
        output_path = Path("exact/random/phi1.txt")
        complete = make_job(
            "exact-random|phi1|4|1|0",
            stage="exact-random",
            output_path=output_path,
            timeout_seconds=30,
        )
        timed_out = make_job(
            "exact-random|phi1|4|1|1",
            stage="exact-random",
            output_path=output_path,
            timeout_seconds=30,
        )
        row = "4 1 - 0 - 0.1 0.2 1 0 1"

        with TemporaryDirectory() as temporary:
            output_directory = Path(temporary)
            state_path = output_directory / runner.STATE_FILENAME
            connection = runner.connect_checkpoint(state_path)
            try:
                runner.initialize_checkpoint(
                    connection, "fingerprint", "source", "config"
                )
                runner.record_attempt(
                    connection,
                    complete,
                    runner.Attempt("complete", row, 0.3, ""),
                )
                runner.record_attempt(
                    connection,
                    timed_out,
                    runner.Attempt("timeout", None, 30.0, "Z3 timeout"),
                )
                records = runner.load_records(connection)
            finally:
                connection.close()

            runner.materialize_results(
                [complete, timed_out], records, output_directory
            )
            self.assertEqual(
                (output_directory / output_path).read_text(encoding="utf-8"),
                row + "\n",
            )
            report = (output_directory / "timeouts.tsv").read_text(
                encoding="utf-8"
            )
            self.assertIn(timed_out.job_id, report)
            self.assertIn("\t30\tZ3 timeout\n", report)

    def test_filtered_timeout_retry_survives_interruption(self) -> None:
        timeout_job = make_job(
            "exact-random|phi2|32|8|37",
            stage="exact-random",
            output_path=Path("exact/random/phi2.txt"),
            timeout_seconds=120,
        )
        selected_job = make_job(
            "exact-random|phi5|4|1|0|grid-2d-v1",
            stage="exact-random",
            output_path=Path("exact/random/phi5.txt"),
            timeout_seconds=120,
        )
        excluded_job = make_job(
            "exact-random|phi5|16|1|0|grid-2d-v1",
            stage="exact-random",
            output_path=Path("exact/random/phi5.txt"),
            timeout_seconds=120,
        )
        jobs = [timeout_job, selected_job, excluded_job]
        rows = (
            "32 8 - 37 - 0.1 0.2 1 0 1",
            "4 1 - 0 - 0.1 0.2 1 0 1",
        )

        with TemporaryDirectory() as temporary:
            output_directory = Path(temporary)
            args = runner.parse_args(
                [
                    "--only",
                    "exact-random",
                    "--formula",
                    "phi5",
                    "--duration",
                    "4",
                    "--rerun-timeouts",
                    "--output-dir",
                    str(output_directory),
                ]
            )
            state_path = output_directory / runner.STATE_FILENAME
            connection = runner.connect_checkpoint(state_path)
            try:
                runner.initialize_checkpoint(
                    connection, "fingerprint", "source", "config"
                )
                runner.record_attempt(
                    connection,
                    timeout_job,
                    runner.Attempt("timeout", None, 120.0, "Z3 timeout"),
                )
            finally:
                connection.close()

            with mock.patch.object(
                runner, "build_plan", return_value=jobs
            ), mock.patch.object(
                runner,
                "experiment_fingerprint",
                return_value=("fingerprint", "source", "config"),
            ), mock.patch.object(
                runner.importlib.util, "find_spec", return_value=object()
            ), mock.patch.object(
                runner, "execute_job", side_effect=KeyboardInterrupt()
            ):
                with contextlib.redirect_stdout(io.StringIO()):
                    self.assertEqual(runner.run(args), 130)

            connection = runner.connect_checkpoint(state_path)
            try:
                records = runner.load_records(connection)
            finally:
                connection.close()
            self.assertEqual(records[timeout_job.job_id]["status"], "timeout")

            with mock.patch.object(
                runner, "build_plan", return_value=jobs
            ), mock.patch.object(
                runner,
                "experiment_fingerprint",
                return_value=("fingerprint", "source", "config"),
            ), mock.patch.object(
                runner.importlib.util, "find_spec", return_value=object()
            ), mock.patch.object(
                runner,
                "execute_job",
                side_effect=[
                    runner.Attempt("complete", rows[0], 0.3, ""),
                    runner.Attempt("complete", rows[1], 0.3, ""),
                ],
            ) as execute:
                with contextlib.redirect_stdout(io.StringIO()):
                    self.assertEqual(runner.run(args), 0)
                self.assertEqual(execute.call_count, 2)

            connection = runner.connect_checkpoint(state_path)
            try:
                records = runner.load_records(connection)
            finally:
                connection.close()
            self.assertEqual(records[timeout_job.job_id]["status"], "complete")
            self.assertEqual(records[selected_job.job_id]["status"], "complete")
            self.assertNotIn(excluded_job.job_id, records)

    def test_checkpoint_rejects_changed_configuration(self) -> None:
        with TemporaryDirectory() as temporary:
            state_path = Path(temporary) / runner.STATE_FILENAME
            connection = runner.connect_checkpoint(state_path)
            runner.initialize_checkpoint(
                connection, "first", "source", '{"exact_repetitions":3}'
            )
            connection.close()

            connection = runner.connect_checkpoint(state_path)
            try:
                with self.assertRaisesRegex(
                    RuntimeError, "different sources or settings"
                ):
                    runner.initialize_checkpoint(
                        connection,
                        "second",
                        "source",
                        '{"exact_repetitions":4}',
                    )
            finally:
                connection.close()

    def test_checkpoint_adopts_timeout_increase_for_retry(self) -> None:
        with TemporaryDirectory() as temporary:
            state_path = Path(temporary) / runner.STATE_FILENAME
            old_config = self.timeout_config(120)
            new_config = self.timeout_config(240)
            connection = runner.connect_checkpoint(state_path)
            runner.initialize_checkpoint(
                connection, "first", "source", old_config
            )
            try:
                changed = runner.initialize_offline_checkpoint(
                    connection,
                    "second",
                    "source",
                    new_config,
                    allow_timeout_increase=True,
                )
                metadata = dict(
                    connection.execute("SELECT key, value FROM metadata")
                )
            finally:
                connection.close()

            self.assertEqual(changed, (False, True))
            self.assertEqual(metadata["fingerprint"], "second")
            self.assertEqual(metadata["config"], new_config)

    def test_checkpoint_adopts_case_study_repetition_change(self) -> None:
        with TemporaryDirectory() as temporary:
            state_path = Path(temporary) / runner.STATE_FILENAME
            previous_config = "{\"case_study_repetitions\":3}"
            current_config = "{\"case_study_repetitions\":10}"
            connection = runner.connect_checkpoint(state_path)
            runner.initialize_checkpoint(
                connection, "first", "source", previous_config
            )
            try:
                changed = runner.initialize_offline_checkpoint(
                    connection,
                    "second",
                    "source",
                    current_config,
                )
                metadata = dict(
                    connection.execute("SELECT key, value FROM metadata")
                )
            finally:
                connection.close()

            self.assertEqual(changed, (True, False))
            self.assertEqual(metadata["fingerprint"], "second")
            self.assertEqual(metadata["config"], current_config)

    def test_timeout_migration_rejects_other_setting_changes(self) -> None:
        self.assertFalse(
            runner.compatible_timeout_increase(
                self.timeout_config(120, exact_repetitions=3),
                self.timeout_config(240, exact_repetitions=4),
            )
        )
        self.assertFalse(
            runner.compatible_timeout_increase(
                self.timeout_config(240),
                self.timeout_config(120),
            )
        )


class ExecutionTests(unittest.TestCase):
    def test_build_selects_release_configuration(self) -> None:
        with mock.patch.object(runner.subprocess, "run") as run_process:
            runner.configure_and_build(
                ("approximate-random",),
                Path("/tmp/sttt-offline-build"),
                2,
                False,
            )
        self.assertEqual(run_process.call_count, 2)
        build_command = run_process.call_args_list[1].args[0]
        config_index = build_command.index("--config")
        self.assertEqual(build_command[config_index + 1], "Release")

    def test_main_reports_sqlite_errors(self) -> None:
        with mock.patch.object(
            runner,
            "run",
            side_effect=runner.sqlite3.OperationalError("broken checkpoint"),
        ), contextlib.redirect_stderr(io.StringIO()) as errors:
            self.assertEqual(runner.main([]), 1)
        self.assertIn("error: broken checkpoint", errors.getvalue())

    def test_exit_124_is_only_an_smt_timeout_and_has_no_wall_timeout(self) -> None:
        exact = make_job(
            "exact-random|phi1|4|1|0",
            stage="exact-random",
            timeout_seconds=30,
        )
        approximate = make_job("approximate-random|phi1|4|1|0")
        process = runner.subprocess.CompletedProcess(
            args=["offline-monitor"],
            returncode=124,
            stdout="",
            stderr="",
        )

        with mock.patch.object(
            runner, "command_for_job", return_value=["offline-monitor"]
        ), mock.patch.object(
            runner.subprocess, "run", return_value=process
        ) as run_process:
            exact_attempt = runner.execute_job(exact, Path("unused-build"))
            approximate_attempt = runner.execute_job(
                approximate, Path("unused-build")
            )

        self.assertEqual(exact_attempt.status, "timeout")
        self.assertEqual(approximate_attempt.status, "error")
        self.assertEqual(run_process.call_count, 2)
        for call in run_process.call_args_list:
            self.assertNotIn("timeout", call.kwargs)


class ResultValidationTests(unittest.TestCase):
    def test_rejects_malformed_result_rows(self) -> None:
        job = make_job("approximate-random|phi1|4|1|0")
        malformed_rows = {
            "multiple rows": "4 1 0 0 7 0.1 1\n4 1 0 0 7 0.1 1\n",
            "wrong column count": "4 1 0 0 7 0.1\n",
            "wrong identity": "8 1 0 0 7 0.1 1\n",
            "non-numeric timing": "4 1 0 0 7 slow 1\n",
            "non-finite timing": "4 1 0 0 7 nan 1\n",
            "negative timing": "4 1 0 0 7 -0.1 1\n",
        }
        with TemporaryDirectory() as temporary:
            result_path = Path(temporary) / job.temporary_path
            for description, contents in malformed_rows.items():
                with self.subTest(description=description):
                    result_path.write_text(contents, encoding="utf-8")
                    with self.assertRaises(RuntimeError):
                        runner.validate_result(job, result_path)

    def test_accepts_and_normalizes_one_valid_row(self) -> None:
        job = make_job("approximate-random|phi1|4|1|0")
        with TemporaryDirectory() as temporary:
            result_path = Path(temporary) / job.temporary_path
            result_path.write_text("  4  1  0  0  7  0.1  1  \n", encoding="utf-8")
            self.assertEqual(
                runner.validate_result(job, result_path),
                "4 1 0 0 7 0.1 1",
            )

    def test_completed_timing_uses_exact_maximum_result_value(self) -> None:
        cases = (
            (
                make_job("approximate-random|phi1|4|1|0"),
                "4 1 0 0 7 0.125000 1",
                "avg 0.125000s",
            ),
            (
                make_job(
                    "approximate-mutual-separation|adm|2|1",
                    stage="approximate-mutual-separation",
                ),
                "1 2 1 adm 0.25 7 1 1",
                "avg 0.25s",
            ),
            (
                make_job(
                    "exact-random|phi1|4|1|0", stage="exact-random"
                ),
                "4 1 - 0 - 0.1 0.2 1 0 1",
                "avg 0.2s",
            ),
            (
                make_job(
                    "exact-mutual-separation|2|1",
                    stage="exact-mutual-separation",
                ),
                "1 2 1 0.4 1 0.5 0 1",
                "avg 0.5s",
            ),
        )
        for job, row, expected in cases:
            with self.subTest(stage=job.stage):
                self.assertEqual(
                    runner.format_attempt_timing(
                        "complete",
                        row,
                        99.0,
                        runner.result_timing_columns(job),
                        "s",
                    ),
                    expected,
                )

    def test_incomplete_timing_is_labeled_as_wall_time(self) -> None:
        self.assertEqual(
            runner.format_attempt_timing(
                "timeout", None, 30.125, (5, 6), "s"
            ),
            "wall 30.125s",
        )


if __name__ == "__main__":
    unittest.main()
