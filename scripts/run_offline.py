#!/usr/bin/env python3
"""Run every offline experiment with resumable configuration checkpoints."""

from __future__ import annotations

import argparse
import fcntl
import importlib.util
import json
import math
import sqlite3
import subprocess
import sys
import tempfile
import time
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

SCRIPT_DIRECTORY = Path(__file__).resolve().parent
if str(SCRIPT_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIRECTORY))

from resumable_runner import (  # noqa: E402
    acquire_lock,
    atomic_write,
    build_context,
    connect_checkpoint,
    diagnostic_text,
    digest_files,
    ensure_fresh_outputs,
    experiment_fingerprint as shared_experiment_fingerprint,
    format_attempt_timing,
    initialize_checkpoint,
    load_records,
    materialize_rows,
    positive_integer,
    record_attempt,
    reset_outputs,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BUILD_DIRECTORY = PROJECT_ROOT / "build"
DEFAULT_OUTPUT_DIRECTORY = PROJECT_ROOT / "results/offline"
STATE_FILENAME = "progress.sqlite3"
LOCK_FILENAME = ".run_offline.lock"
PLAN_VERSION = 1
EXACT_TIMED_RANDOM_BACKEND = "grid-2d-v1"
EXACT_TIMED_RANDOM_FORMULAS = {"phi5", "phi6"}
EXACT_CASE_STUDY_BACKEND = "rv-full-flow-v1"

STAGES = (
    "approximate-random",
    "exact-random",
    "approximate-mutual-separation",
    "approximate-water-tanks",
    "exact-mutual-separation",
    "exact-water-tanks",
)
STAGE_DESCRIPTIONS = {
    "approximate-random": "approximate monitor on random traces",
    "approximate-mutual-separation": "approximate mutual-separation case study",
    "approximate-water-tanks": "approximate water-tank case study",
    "exact-random": "exact SMT monitor on random traces",
    "exact-mutual-separation": "exact mutual-separation case study",
    "exact-water-tanks": "exact water-tank case study",
}
APPROXIMATE_TARGETS = {
    "approximate-random": "sttt_offline_approximate_random",
    "approximate-mutual-separation":
        "sttt_offline_approximate_mutual_separation",
    "approximate-water-tanks": "sttt_offline_approximate_water_tanks",
}
APPROXIMATE_BINARIES = {
    stage: target for stage, target in APPROXIMATE_TARGETS.items()
}
VARIANT_SUFFIXES = {
    "adm": "",
    "adm-f": "_adm_f",
    "adm-fr": "_adm_fr",
    "adm-c": "_adm_c",
    "adm-cr": "_adm_cr",
}


@dataclass(frozen=True)
class Job:
    job_id: str
    stage: str
    label: str
    entrypoint: str
    arguments: tuple[str, ...]
    output_path: Path
    temporary_path: Path
    timeout_seconds: int | None


@dataclass(frozen=True)
class Attempt:
    status: str
    row: str | None
    elapsed_seconds: float
    message: str


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--only",
        action="append",
        choices=STAGES,
        metavar="STAGE",
        help="run only this stage; may be supplied more than once",
    )
    parser.add_argument(
        "--formula",
        action="append",
        choices=tuple(f"phi{index}" for index in range(1, 7)),
        dest="formulas",
        help="random formula to run; may be supplied more than once",
    )
    parser.add_argument(
        "--duration",
        action="append",
        choices=(4, 8, 16, 32),
        dest="durations",
        type=int,
        help="random duration to run; may be supplied more than once",
    )
    parser.add_argument(
        "--approximate-repetitions",
        type=positive_integer,
        default=100,
        help="timed repetitions for approximate random monitors (default: 100)",
    )
    parser.add_argument(
        "--exact-repetitions",
        type=positive_integer,
        default=3,
        help="timed repetitions for exact random monitors (default: 3)",
    )
    parser.add_argument(
        "--case-study-repetitions",
        type=positive_integer,
        default=10,
        help="timed repetitions for all case-study monitors (default: 10)",
    )
    parser.add_argument(
        "--random-timeout",
        type=positive_integer,
        default=600,
        help="SMT timeout per random solver check (default: 600)",
    )
    parser.add_argument(
        "--mutual-separation-timeout",
        type=positive_integer,
        default=360,
        help="SMT timeout per mutual-separation solver check (default: 360)",
    )
    parser.add_argument(
        "--water-tanks-timeout",
        type=positive_integer,
        default=120,
        help="SMT timeout per water-tank solver check (default: 120)",
    )
    parser.add_argument(
        "--build-dir",
        type=Path,
        default=DEFAULT_BUILD_DIRECTORY,
        help="CMake build directory",
    )
    parser.add_argument(
        "--build-jobs",
        type=positive_integer,
        default=2,
        help="parallel CMake build jobs (default: 2)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIRECTORY,
        help="generated result and checkpoint directory",
    )
    parser.add_argument(
        "--no-build",
        action="store_true",
        help="use existing approximate binaries without invoking CMake",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="show the selected plan without building or running it",
    )
    parser.add_argument(
        "--status",
        action="store_true",
        help="show checkpoint progress without changing it",
    )
    parser.add_argument(
        "--restart",
        action="store_true",
        help="discard the checkpoint and start all experiments again",
    )
    parser.add_argument(
        "--rerun-timeouts",
        action="store_true",
        help=(
            "retry all timed-out jobs in selected stages, including jobs "
            "outside formula and duration filters"
        ),
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="list offline stages and exit",
    )
    args = parser.parse_args(argv)

    if args.status and (args.restart or args.rerun_timeouts):
        parser.error("--status cannot be combined with --restart or --rerun-timeouts")
    if args.dry_run and (args.restart or args.rerun_timeouts):
        parser.error("--dry-run cannot be combined with --restart or --rerun-timeouts")

    args.build_dir = args.build_dir.resolve()
    args.output_dir = args.output_dir.resolve()
    return args


def selected_stages(args: argparse.Namespace) -> tuple[str, ...]:
    if not args.only:
        return STAGES
    selected = set(args.only)
    return tuple(stage for stage in STAGES if stage in selected)

def jobs_for_stages(
    plan: Sequence[Job], stages: Sequence[str]
) -> list[Job]:
    """Return jobs grouped in the requested canonical stage order."""
    return [job for stage in stages for job in plan if job.stage == stage]


def select_jobs(
    plan: Sequence[Job],
    stages: Sequence[str],
    args: argparse.Namespace,
    additional_job_ids: set[str] | None = None,
) -> list[Job]:
    """Apply random formula/duration filters in canonical stage order."""
    formulas = set(args.formulas or ())
    durations = set(args.durations or ())
    additional = additional_job_ids or set()
    selected: list[Job] = []

    for job in jobs_for_stages(plan, stages):
        if job.job_id in additional or not job.stage.endswith("-random"):
            selected.append(job)
            continue

        parts = job.job_id.split("|")
        formula = parts[1]
        duration = int(parts[2])
        if formulas and formula not in formulas:
            continue
        if durations and duration not in durations:
            continue
        selected.append(job)

    return selected


def print_stages() -> None:
    for stage in STAGES:
        print(f"{stage}: {STAGE_DESCRIPTIONS[stage]}")


def build_plan(args: argparse.Namespace) -> list[Job]:
    jobs: list[Job] = []

    def add_job(
        stage: str,
        job_id: str,
        label: str,
        entrypoint: str,
        arguments: tuple[str, ...],
        output_path: Path,
        temporary_path: Path,
        repetitions: int,
        timeout_seconds: int | None,
    ) -> None:
        timeout = timeout_seconds if stage.startswith("exact-") else None
        extra_arguments = ["--repeat", str(repetitions)]
        if timeout is not None:
            extra_arguments.extend(("--timeout", str(timeout)))
        jobs.append(
            Job(
                job_id=job_id,
                stage=stage,
                label=label,
                entrypoint=entrypoint,
                arguments=(*arguments, *extra_arguments),
                output_path=output_path,
                temporary_path=temporary_path,
                timeout_seconds=timeout,
            )
        )

    random_binary = APPROXIMATE_BINARIES["approximate-random"]
    for formula in (f"phi{index}" for index in range(1, 7)):
        for duration in (4, 8, 16, 32):
            for epsilon in (1, 2, 4, 8):
                if epsilon > duration:
                    continue
                for sample in range(100):
                    label = (
                        f"{formula} duration={duration} epsilon={epsilon} "
                        f"sample={sample}"
                    )
                    job_id = (
                        f"approximate-random|{formula}|{duration}|{epsilon}|"
                        f"{sample}"
                    )
                    add_job(
                        "approximate-random",
                        job_id,
                        label,
                        random_binary,
                        (
                            "--formula",
                            formula,
                            "--duration",
                            str(duration),
                            "--epsilon",
                            str(epsilon),
                            "--sample",
                            str(sample),
                            "--data-dir",
                            str(PROJECT_ROOT / "data/signals"),
                        ),
                        Path("approximate/random") / f"{formula}.txt",
                        Path(f"{formula}.txt"),
                        args.approximate_repetitions,
                        None,
                    )

    case_definitions = (
        (
            "approximate-mutual-separation",
            "mutual_separation",
            ("adm", "adm-f", "adm-fr"),
            (1, 2, 3, 4, 5),
            "data/case_studies/uav",
        ),
        (
            "approximate-water-tanks",
            "water_tanks",
            ("adm", "adm-f", "adm-fr", "adm-c", "adm-cr"),
            (1, 2, 4, 8),
            "data/case_studies/water_tanks",
        ),
    )
    for stage, subject, variants, epsilons, data_path in case_definitions:
        for variant in variants:
            filename = f"{subject}{VARIANT_SUFFIXES[variant]}.txt"
            for agents in (2, 3, 4):
                for epsilon in epsilons:
                    label = (
                        f"{variant} agents={agents} epsilon={epsilon}"
                    )
                    job_id = (
                        f"{stage}|{variant}|{agents}|{epsilon}"
                    )
                    add_job(
                        stage,
                        job_id,
                        label,
                        APPROXIMATE_BINARIES[stage],
                        (
                            "--agents",
                            str(agents),
                            "--epsilon-multiplier",
                            str(epsilon),
                            "--variant",
                            variant,
                            "--data-dir",
                            str(PROJECT_ROOT / data_path),
                        ),
                        Path("approximate/case_studies") / filename,
                        Path(filename),
                        args.case_study_repetitions,
                        None,
                    )

    exact_entrypoint = str(PROJECT_ROOT / "benchmarks/offline/exact/run.py")
    for formula in (f"phi{index}" for index in range(1, 7)):
        timed_formula = formula in EXACT_TIMED_RANDOM_FORMULAS
        durations = (4, 8) if timed_formula else (4, 8, 16, 32)
        for duration in durations:
            for epsilon in (1, 2, 4, 8):
                if epsilon > duration:
                    continue
                for sample in range(100):
                    label = (
                        f"{formula} duration={duration} epsilon={epsilon} "
                        f"sample={sample}"
                    )
                    job_id = (
                        f"exact-random|{formula}|{duration}|{epsilon}|{sample}"
                    )
                    if timed_formula:
                        job_id += f"|{EXACT_TIMED_RANDOM_BACKEND}"
                    add_job(
                        "exact-random",
                        job_id,
                        label,
                        exact_entrypoint,
                        (
                            "random",
                            "--formula",
                            formula,
                            "--duration",
                            str(duration),
                            "--epsilon",
                            str(epsilon),
                            "--sample",
                            str(sample),
                        ),
                        Path("exact/random") / f"{formula}.txt",
                        Path("random") / f"{formula}.txt",
                        args.exact_repetitions,
                        args.random_timeout,
                    )

    exact_cases = (
        (
            "exact-mutual-separation",
            "mutual-separation",
            "mutual_separation.txt",
            (1, 2, 3, 4, 5),
            args.mutual_separation_timeout,
        ),
        (
            "exact-water-tanks",
            "water-tanks",
            "water_tanks.txt",
            (1, 2, 4, 8),
            args.water_tanks_timeout,
        ),
    )
    for stage, subject, filename, epsilons, timeout in exact_cases:
        for agents in (2, 3, 4):
            for epsilon in epsilons:
                label = f"agents={agents} epsilon={epsilon}"
                job_id = (
                    f"{stage}|{agents}|{epsilon}|{EXACT_CASE_STUDY_BACKEND}"
                )
                add_job(
                    stage,
                    job_id,
                    label,
                    exact_entrypoint,
                    (
                        subject,
                        "--agents",
                        str(agents),
                        "--epsilon",
                        str(epsilon),
                    ),
                    Path("exact/case_studies") / filename,
                    Path("case_studies") / filename,
                    args.case_study_repetitions,
                    timeout,
                )

    return jobs


def relevant_files() -> list[Path]:
    paths = {
        Path(__file__).resolve(),
        SCRIPT_DIRECTORY / "resumable_runner.py",
        PROJECT_ROOT / "CMakeLists.txt",
        PROJECT_ROOT / "requirements-smt.txt",
    }
    for directory in (
        PROJECT_ROOT / "include",
        PROJECT_ROOT / "benchmarks/offline",
        PROJECT_ROOT / "data",
    ):
        if not directory.exists():
            continue
        for path in directory.rglob("*"):
            if (
                path.is_file()
                and "__pycache__" not in path.parts
                and path.suffix not in (".pyc", ".pyo", ".orig", ".rej")
                and not path.name.endswith("~")
            ):
                paths.add(path.resolve())
    return sorted(paths, key=lambda path: str(path.relative_to(PROJECT_ROOT)))


def source_digest() -> str:
    return digest_files(PROJECT_ROOT, relevant_files())


def experiment_config(args: argparse.Namespace) -> dict[str, object]:
    return {
        "plan_version": PLAN_VERSION,
        "build": build_context(args.build_dir),
        "no_build": args.no_build,
        "approximate_repetitions": args.approximate_repetitions,
        "exact_repetitions": args.exact_repetitions,
        "case_study_repetitions": args.case_study_repetitions,
        "random_timeout": args.random_timeout,
        "mutual_separation_timeout": args.mutual_separation_timeout,
        "water_tanks_timeout": args.water_tanks_timeout,
    }


def experiment_fingerprint(
    args: argparse.Namespace,
) -> tuple[str, str, str]:
    return shared_experiment_fingerprint(
        PROJECT_ROOT,
        relevant_files(),
        experiment_config(args),
    )


def output_paths(plan: Sequence[Job], output_directory: Path) -> set[Path]:
    paths = {output_directory / job.output_path for job in plan}
    paths.add(output_directory / "timeouts.tsv")
    return paths


def command_for_job(job: Job, build_directory: Path, output: Path) -> list[str]:
    if job.stage.startswith("approximate-"):
        executable = build_directory / "bin" / job.entrypoint
        command = [str(executable)]
    else:
        command = [sys.executable, job.entrypoint]
    return [
        *command,
        *job.arguments,
        "--output-dir",
        str(output),
        "--quiet",
    ]


def result_timing_columns(job: Job) -> tuple[int, ...]:
    if job.stage == "approximate-random":
        return (5,)
    if job.stage.startswith("approximate-"):
        return (4,)
    if job.stage == "exact-random":
        return (5, 6)
    return (3, 5)


def validate_result(job: Job, result_path: Path) -> str:
    if not result_path.is_file():
        raise RuntimeError(f"runner did not create {job.temporary_path}")

    rows = [
        line.strip()
        for line in result_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(rows) != 1:
        raise RuntimeError(f"expected one result row, found {len(rows)}")

    fields = rows[0].split()
    parts = job.job_id.split("|")
    if job.stage == "approximate-random":
        expected_columns = 7
        expected_values = (
            (0, parts[2]),
            (1, parts[3]),
            (2, "0"),
            (3, parts[4]),
        )
    elif job.stage.startswith("approximate-"):
        expected_columns = 8
        expected_values = (
            (0, "1"),
            (1, parts[2]),
            (2, parts[3]),
        )
    elif job.stage == "exact-random":
        expected_columns = 10
        expected_values = (
            (0, parts[2]),
            (1, parts[3]),
            (2, "-"),
            (3, parts[4]),
            (4, "-"),
        )
    else:
        expected_columns = 8
        expected_values = (
            (0, "1"),
            (1, parts[1]),
            (2, parts[2]),
        )

    if len(fields) != expected_columns:
        raise RuntimeError(
            f"expected {expected_columns} result fields, found {len(fields)}"
        )
    for index, expected in expected_values:
        if fields[index] != expected:
            raise RuntimeError(
                f"result field {index} is {fields[index]!r}; "
                f"expected {expected!r}"
            )
    for index in result_timing_columns(job):
        try:
            timing = float(fields[index])
        except ValueError as error:
            raise RuntimeError("result contains a non-numeric timing") from error
        if not math.isfinite(timing) or timing < 0:
            raise RuntimeError("result contains an invalid timing")
    return " ".join(fields)


def execute_job(job: Job, build_directory: Path) -> Attempt:
    started = time.perf_counter()
    with tempfile.TemporaryDirectory(prefix="sttt-offline-") as temporary:
        output = Path(temporary)
        command = command_for_job(job, build_directory, output)
        process = subprocess.run(
            command,
            cwd=PROJECT_ROOT,
            capture_output=True,
            text=True,
            check=False,
        )

        elapsed = time.perf_counter() - started
        detail = diagnostic_text(process.stderr or process.stdout)
        if process.returncode == 124 and job.timeout_seconds is not None:
            return Attempt(
                "timeout",
                None,
                elapsed,
                detail or f"SMT timeout after {job.timeout_seconds}s",
            )
        if process.returncode != 0:
            return Attempt(
                "error",
                None,
                elapsed,
                detail or f"process exited with status {process.returncode}",
            )

        try:
            row = validate_result(job, output / job.temporary_path)
        except RuntimeError as error:
            return Attempt("error", None, elapsed, str(error))
        return Attempt("complete", row, elapsed, "")


def materialize_results(
    plan: Sequence[Job],
    records: dict[str, sqlite3.Row],
    output_directory: Path,
) -> None:
    materialize_rows(plan, records, output_directory)

    timeout_lines = ["stage\tinstance\ttimeout_seconds\treason\n"]
    for job in plan:
        record = records.get(job.job_id)
        if record is None or record["status"] != "timeout":
            continue
        reason = record["message"].replace("\t", " ").replace("\n", " ")
        timeout = "" if job.timeout_seconds is None else str(job.timeout_seconds)
        timeout_lines.append(
            f"{job.stage}\t{job.job_id}\t{timeout}\t{reason}\n"
        )
    atomic_write(output_directory / "timeouts.tsv", "".join(timeout_lines))


def configure_and_build(
    stages: Sequence[str],
    build_directory: Path,
    build_jobs: int,
    no_build: bool,
) -> None:
    approximate = [
        stage for stage in stages if stage in APPROXIMATE_TARGETS
    ]
    if not approximate:
        return

    binaries = [
        build_directory / "bin" / APPROXIMATE_BINARIES[stage]
        for stage in approximate
    ]
    if no_build:
        missing = [str(path) for path in binaries if not path.is_file()]
        if missing:
            raise RuntimeError(
                "missing approximate binaries: "
                + ", ".join(missing)
                + "; omit --no-build to build them"
            )
        return

    subprocess.run(
        [
            "cmake",
            "-S",
            str(PROJECT_ROOT),
            "-B",
            str(build_directory),
            "-DCMAKE_BUILD_TYPE=Release",
        ],
        cwd=PROJECT_ROOT,
        check=True,
    )
    subprocess.run(
        [
            "cmake",
            "--build",
            str(build_directory),
            "--config",
            "Release",
            "--parallel",
            str(build_jobs),
            "--target",
            *(APPROXIMATE_TARGETS[stage] for stage in approximate),
        ],
        cwd=PROJECT_ROOT,
        check=True,
    )


def print_plan(
    args: argparse.Namespace,
    plan: Sequence[Job],
    stages: Sequence[str],
) -> None:
    selected = select_jobs(plan, stages, args)
    print("Offline experiment plan:")
    counts = Counter(job.stage for job in selected)
    for stage in stages:
        print(f"  {stage}: {counts[stage]} instances")
    print(f"Total: {len(selected)} instances")
    print(
        "Repetitions per instance: "
        f"random approximate={args.approximate_repetitions}, "
        f"random exact={args.exact_repetitions}, "
        f"case studies={args.case_study_repetitions}"
    )
    print(
        "SMT timeouts per solver check: "
        f"random={args.random_timeout}s, "
        f"mutual separation={args.mutual_separation_timeout}s, "
        f"water tanks={args.water_tanks_timeout}s"
    )
    print("Approximate monitors have no timeout.")


def initialize_offline_checkpoint(
    connection: sqlite3.Connection,
    fingerprint: str,
    source: str,
    config: str,
    allow_timeout_increase: bool = False,
) -> tuple[bool, bool]:
    """Adopt source or setting changes while preserving checkpoint rows."""
    metadata = dict(connection.execute("SELECT key, value FROM metadata"))
    previous = metadata.get("fingerprint")
    if previous is None:
        initialize_checkpoint(connection, fingerprint, source, config)
        return False, False
    if previous == fingerprint:
        return False, False
    previous_config = metadata.get("config")
    timeout_increased = (
        allow_timeout_increase
        and previous_config is not None
        and compatible_timeout_increase(previous_config, config)
    )
    sources_changed = metadata.get("source_digest") != source
    checkpoint_changed = sources_changed or (
        previous_config != config and not timeout_increased
    )
    with connection:
        connection.executemany(
            "UPDATE metadata SET value = ? WHERE key = ?",
            (
                (fingerprint, "fingerprint"),
                (source, "source_digest"),
                (config, "config"),
            ),
        )
    return checkpoint_changed, timeout_increased

def compatible_timeout_increase(previous: str, current: str) -> bool:
    """Return whether configs differ only by one or more timeout increases."""
    try:
        previous_config = json.loads(previous)
        current_config = json.loads(current)
    except (json.JSONDecodeError, TypeError):
        return False
    if not isinstance(previous_config, dict) or not isinstance(
        current_config, dict
    ):
        return False
    if previous_config.keys() != current_config.keys():
        return False

    timeout_keys = {
        "random_timeout",
        "mutual_separation_timeout",
        "water_tanks_timeout",
    }
    increased = False
    for key, previous_value in previous_config.items():
        current_value = current_config[key]
        if key not in timeout_keys:
            if previous_value != current_value:
                return False
            continue
        if not isinstance(previous_value, int) or not isinstance(
            current_value, int
        ):
            return False
        if current_value < previous_value:
            return False
        increased = increased or current_value > previous_value
    return increased


def status_report(
    state_path: Path,
    plan: Sequence[Job],
    stages: Sequence[str],
    args: argparse.Namespace,
) -> int:
    if not state_path.is_file():
        print(f"No checkpoint exists at {state_path}")
        return 0

    try:
        connection = sqlite3.connect(
            f"file:{state_path}?mode=ro",
            uri=True,
        )
        connection.row_factory = sqlite3.Row
        metadata = dict(connection.execute("SELECT key, value FROM metadata"))
        records = load_records(connection)
    except sqlite3.DatabaseError as error:
        raise RuntimeError(f"cannot read checkpoint {state_path}: {error}") from error
    finally:
        if "connection" in locals():
            connection.close()

    fingerprint, _, _ = experiment_fingerprint(args)
    if metadata.get("fingerprint") != fingerprint:
        print("Warning: checkpoint sources or settings differ from this invocation.")

    selected = select_jobs(plan, stages, args)
    print(f"Checkpoint: {state_path}")
    for stage in stages:
        stage_jobs = [job for job in selected if job.stage == stage]
        counts = Counter(
            records[job.job_id]["status"]
            for job in stage_jobs
            if job.job_id in records
        )
        not_started = len(stage_jobs) - sum(counts.values())
        print(
            f"{stage}: {counts['complete']} complete, "
            f"{counts['timeout']} timed out, {counts['error']} failed, "
            f"{not_started} not started"
        )
    return 0


def print_summary(
    jobs: Sequence[Job],
    records: dict[str, sqlite3.Row],
    output_directory: Path,
) -> None:
    counts = Counter(
        records[job.job_id]["status"]
        for job in jobs
        if job.job_id in records
    )
    not_started = len(jobs) - sum(counts.values())
    print(
        f"Summary: {counts['complete']} complete, "
        f"{counts['timeout']} timed out, {counts['error']} failed, "
        f"{not_started} not started."
    )
    print(f"Results: {output_directory}")


def _run(args: argparse.Namespace) -> int:
    stages = selected_stages(args)
    plan = build_plan(args)

    if args.list:
        print_stages()
        return 0
    if args.dry_run:
        print_plan(args, plan, stages)
        return 0

    state_path = args.output_dir / STATE_FILENAME
    if args.status:
        return status_report(state_path, plan, stages, args)

    paths = output_paths(plan, args.output_dir)

    if args.restart:
        reset_outputs(state_path, paths)
    else:
        ensure_fresh_outputs(state_path, paths)

    fingerprint, source, config = experiment_fingerprint(args)
    connection = connect_checkpoint(state_path)
    try:
        checkpoint_changed, timeout_increased = initialize_offline_checkpoint(
            connection,
            fingerprint,
            source,
            config,
            allow_timeout_increase=args.rerun_timeouts,
        )
        if checkpoint_changed:
            print(
                "Warning: checkpoint sources or settings changed; preserving "
                "completed instances and continuing."
            )
        if timeout_increased:
            print(
                "Checkpoint SMT timeout increased; preserving completed "
                "instances and retrying timed-out instances."
            )

        records = load_records(connection)
        stage_jobs = jobs_for_stages(plan, stages)
        retry_job_ids = {
            job.job_id
            for job in stage_jobs
            if (
                job.job_id in records
                and records[job.job_id]["status"] in ("timeout", "error")
            )
        }
        timed_out_job_ids = {
            job_id
            for job_id in retry_job_ids
            if records[job_id]["status"] == "timeout"
        }
        selected = select_jobs(
            plan,
            stages,
            args,
            retry_job_ids if args.rerun_timeouts else None,
        )
        materialize_results(plan, records, args.output_dir)

        pending = [
            job
            for job in selected
            if (
                job.job_id not in records
                or records[job.job_id]["status"] not in ("complete", "timeout")
                or (
                    args.rerun_timeouts
                    and job.job_id in timed_out_job_ids
                )
            )
        ]
        pending_stages = tuple(
            stage
            for stage in stages
            if any(job.stage == stage for job in pending)
        )
        if (
            any(stage.startswith("exact-") for stage in pending_stages)
            and importlib.util.find_spec("z3") is None
        ):
            raise RuntimeError(
                "Z3 is not installed; install requirements-smt.txt and retry"
            )
        configure_and_build(
            pending_stages,
            args.build_dir,
            args.build_jobs,
            args.no_build,
        )

        terminal = {
            job_id
            for job_id, record in records.items()
            if record["status"] in ("complete", "timeout")
            and not (
                args.rerun_timeouts
                and job_id in timed_out_job_ids
            )
        }
        finished = sum(job.job_id in terminal for job in selected)
        print(
            f"Running {len(selected)} checkpointed runs; "
            f"{finished} already finished."
        )

        previous_stage: str | None = None
        for job in selected:
            record = records.get(job.job_id)
            if (
                record is not None
                and (
                    record["status"] == "complete"
                    or (
                        record["status"] == "timeout"
                        and not args.rerun_timeouts
                    )
                )
            ):
                continue

            if previous_stage is not None and job.stage != previous_stage:
                materialize_results(plan, records, args.output_dir)
            if job.stage != previous_stage:
                print(f"\n{job.stage}")
                previous_stage = job.stage

            try:
                attempt = execute_job(job, args.build_dir)
            except KeyboardInterrupt:
                materialize_results(plan, records, args.output_dir)
                print("\nInterrupted; completed jobs are checkpointed.")
                print_summary(selected, records, args.output_dir)
                return 130

            records[job.job_id] = record_attempt(connection, job, attempt)
            if attempt.status in ("complete", "timeout"):
                finished += 1

            if (
                sys.stdout.isatty()
                or attempt.status != "complete"
                or finished % 100 == 0
                or finished == len(selected)
            ):
                timing = format_attempt_timing(
                    attempt.status,
                    attempt.row,
                    attempt.elapsed_seconds,
                    result_timing_columns(job),
                    "s",
                )
                print(
                    f"[{finished}/{len(selected)}] {job.label}: "
                    f"{attempt.status} ({timing})",
                    flush=True,
                )

            if attempt.status == "error":
                materialize_results(plan, records, args.output_dir)
                print(f"Error in {job.stage} {job.label}: {attempt.message}")
                print("The failed job will be retried on the next invocation.")
                print_summary(selected, records, args.output_dir)
                return 1

        materialize_results(plan, records, args.output_dir)
        print_summary(selected, records, args.output_dir)
        return 0
    finally:
        connection.close()


def run(args: argparse.Namespace) -> int:
    if args.list or args.dry_run or args.status:
        return _run(args)

    lock = acquire_lock(args.output_dir, LOCK_FILENAME, "offline")
    try:
        return _run(args)
    finally:
        fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
        lock.close()


def main(argv: Sequence[str] | None = None) -> int:
    try:
        return run(parse_args(argv))
    except (
        OSError,
        RuntimeError,
        sqlite3.Error,
        subprocess.CalledProcessError,
    ) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("\nInterrupted.", file=sys.stderr)
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
