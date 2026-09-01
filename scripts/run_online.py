#!/usr/bin/env python3
"""Run every online experiment with resumable configuration checkpoints."""

from __future__ import annotations

import argparse
import csv
import fcntl
import io
import math
import sqlite3
import subprocess
import sys
import tempfile
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from statistics import fmean
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
DEFAULT_OUTPUT_DIRECTORY = PROJECT_ROOT / "results/online"
STATE_FILENAME = "progress.sqlite3"
LOCK_FILENAME = ".run_online.lock"
SUMMARY_FILENAME = "summary.csv"
PLAN_VERSION = 1

MODES = ("incremental", "naive")
FORMULAS = (
    "ac",
    "ad",
    "since",
    "response",
    "response01",
    "response02",
)
PAPER_FORMULAS = {
    "ac": "Psi 1",
    "ad": "Psi 2",
    "since": "Psi 3",
    "response": "Psi 4",
    "response01": "Psi 5",
    "response02": "Psi 6",
}
DURATIONS = (64, 128, 256, 512)
EPSILONS = (2, 4, 8)
CHUNK_SIZES = (4, 8, 16, 32)
SAMPLES = range(100)


@dataclass(frozen=True)
class Job:
    job_id: str
    stage: str
    label: str
    target: str
    arguments: tuple[str, ...]
    output_path: Path
    temporary_path: Path
    mode: str
    formula: str
    duration: int
    epsilon: int
    chunk_size: int
    sample: int


@dataclass(frozen=True)
class Attempt:
    status: str
    row: str | None
    elapsed_seconds: float
    message: str


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mode",
        choices=(*MODES, "all"),
        default="all",
        help="implementation to run (default: all)",
    )
    parser.add_argument(
        "--formula",
        action="append",
        choices=FORMULAS,
        dest="formulas",
        help="formula to run; repeat to select several (default: all)",
    )
    parser.add_argument(
        "--repetitions",
        type=positive_integer,
        default=10,
        help="timed repetitions per configuration (default: 10)",
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
        help="use existing online binaries without invoking CMake",
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
        help="discard the checkpoint and start all online experiments again",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="list selected online CMake targets and exit",
    )
    args = parser.parse_args(argv)

    if args.status and args.restart:
        parser.error("--status cannot be combined with --restart")
    if args.dry_run and args.restart:
        parser.error("--dry-run cannot be combined with --restart")

    args.build_dir = args.build_dir.resolve()
    args.output_dir = args.output_dir.resolve()
    return args


def selected_stages(args: argparse.Namespace) -> tuple[str, ...]:
    modes = MODES if args.mode == "all" else (args.mode,)
    formulas = args.formulas or FORMULAS
    selected = set(formulas)
    return tuple(
        f"{mode}/{formula}"
        for mode in modes
        for formula in FORMULAS
        if formula in selected
    )


def output_name(mode: str, formula: str) -> str:
    prefix = "onlineIncr" if mode == "incremental" else "onlineNaive"
    return f"{prefix}_{formula}.txt"


def build_plan(args: argparse.Namespace) -> list[Job]:
    jobs: list[Job] = []
    data_directory = PROJECT_ROOT / "data/signals"
    for mode in MODES:
        for formula in FORMULAS:
            stage = f"{mode}/{formula}"
            target = f"sttt_online_{mode}_{formula}"
            filename = output_name(mode, formula)
            relative_output = Path(mode) / filename
            for duration in DURATIONS:
                for epsilon in EPSILONS:
                    for chunk_size in CHUNK_SIZES:
                        if chunk_size <= epsilon or chunk_size > duration:
                            continue
                        for sample in SAMPLES:
                            job_id = (
                                f"{mode}|{formula}|{duration}|{epsilon}|"
                                f"{chunk_size}|{sample}"
                            )
                            label = (
                                f"{formula} duration={duration} "
                                f"epsilon={epsilon} chunk={chunk_size} "
                                f"sample={sample}"
                            )
                            jobs.append(
                                Job(
                                    job_id=job_id,
                                    stage=stage,
                                    label=label,
                                    target=target,
                                    arguments=(
                                        "--duration",
                                        str(duration),
                                        "--epsilon",
                                        str(epsilon),
                                        "--chunk",
                                        str(chunk_size),
                                        "--sample",
                                        str(sample),
                                        "--repeat",
                                        str(args.repetitions),
                                        "--data-dir",
                                        str(data_directory),
                                    ),
                                    output_path=relative_output,
                                    temporary_path=relative_output,
                                    mode=mode,
                                    formula=formula,
                                    duration=duration,
                                    epsilon=epsilon,
                                    chunk_size=chunk_size,
                                    sample=sample,
                                )
                            )
    return jobs


def relevant_files() -> list[Path]:
    paths = {
        Path(__file__).resolve(),
        SCRIPT_DIRECTORY / "resumable_runner.py",
        PROJECT_ROOT / "CMakeLists.txt",
    }
    for directory in (
        PROJECT_ROOT / "include",
        PROJECT_ROOT / "benchmarks/online",
        PROJECT_ROOT / "data/signals",
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


def experiment_config(args: argparse.Namespace) -> dict[str, object]:
    return {
        "plan_version": PLAN_VERSION,
        "build": build_context(args.build_dir),
        "no_build": args.no_build,
        "repetitions": args.repetitions,
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
    paths.add(output_directory / SUMMARY_FILENAME)
    return paths


def command_for_job(job: Job, build_directory: Path, output: Path) -> list[str]:
    return [
        str(build_directory / "bin" / job.target),
        *job.arguments,
        "--output-dir",
        str(output / job.mode),
        "--quiet",
    ]


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
    if len(fields) != 6:
        raise RuntimeError(f"expected 6 result fields, found {len(fields)}")
    expected = (
        str(job.duration),
        str(job.epsilon),
        str(job.chunk_size),
        str(job.sample),
    )
    if tuple(fields[:4]) != expected:
        raise RuntimeError(
            f"result key {tuple(fields[:4])!r}; expected {expected!r}"
        )
    try:
        timing = float(fields[4])
    except ValueError as error:
        raise RuntimeError("result contains a non-numeric timing") from error
    if not math.isfinite(timing) or timing < 0:
        raise RuntimeError("result contains an invalid timing")

    verdicts = fields[5]
    expected_verdicts = math.ceil(job.duration / job.chunk_size)
    if len(verdicts) != expected_verdicts or any(
        verdict not in "012" for verdict in verdicts
    ):
        raise RuntimeError(
            "result contains an invalid online verdict sequence"
        )
    return " ".join(fields)


def execute_job(job: Job, build_directory: Path) -> Attempt:
    started = time.perf_counter()
    with tempfile.TemporaryDirectory(prefix="sttt-online-") as temporary:
        output = Path(temporary)
        process = subprocess.run(
            command_for_job(job, build_directory, output),
            cwd=PROJECT_ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
        elapsed = time.perf_counter() - started
        detail = diagnostic_text(process.stderr or process.stdout)
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


def parse_completed_row(job: Job, row: str) -> tuple[float, str]:
    fields = row.split()
    if len(fields) != 6:
        raise RuntimeError(f"checkpoint contains an invalid row for {job.job_id}")
    try:
        timing = float(fields[4])
    except ValueError as error:
        raise RuntimeError(
            f"checkpoint contains an invalid timing for {job.job_id}"
        ) from error
    return timing, fields[5]


def materialize_summary(
    plan: Sequence[Job],
    records: dict[str, sqlite3.Row],
    output_directory: Path,
) -> None:
    samples: dict[
        tuple[str, int, int, int],
        dict[str, dict[int, tuple[float, str]]],
    ] = defaultdict(lambda: defaultdict(dict))
    for job in plan:
        if job.formula not in PAPER_FORMULAS:
            continue
        record = records.get(job.job_id)
        if (
            record is None
            or record["status"] != "complete"
            or record["result_row"] is None
        ):
            continue
        key = (job.formula, job.duration, job.epsilon, job.chunk_size)
        samples[key][job.mode][job.sample] = parse_completed_row(
            job, record["result_row"]
        )

    buffer = io.StringIO()
    writer = csv.writer(buffer, lineterminator="\n")
    writer.writerow(
        (
            "formula",
            "d",
            "eps",
            "Delta",
            "speedup",
            "incr wall time",
            "time per portion",
        )
    )
    for formula, label in PAPER_FORMULAS.items():
        for duration in DURATIONS:
            for epsilon in EPSILONS:
                for chunk_size in CHUNK_SIZES:
                    if chunk_size <= epsilon or chunk_size > duration:
                        continue
                    key = (formula, duration, epsilon, chunk_size)
                    modes = samples.get(key, {})
                    incremental = modes.get("incremental", {})
                    naive = modes.get("naive", {})
                    if len(incremental) != len(SAMPLES) or len(naive) != len(SAMPLES):
                        continue

                    ratios = []
                    incremental_times = []
                    for sample in SAMPLES:
                        incremental_time, incremental_verdicts = incremental[sample]
                        naive_time, naive_verdicts = naive[sample]
                        if incremental_verdicts != naive_verdicts:
                            raise RuntimeError(
                                "online implementations disagree for "
                                f"{formula} duration={duration} epsilon={epsilon} "
                                f"chunk={chunk_size} sample={sample}"
                            )
                        if incremental_time <= 0:
                            raise RuntimeError(
                                f"incremental timing is not positive for {key}"
                            )
                        ratios.append(naive_time / incremental_time)
                        incremental_times.append(incremental_time)

                    incremental_mean = fmean(incremental_times)
                    writer.writerow(
                        (
                            label,
                            duration,
                            epsilon,
                            chunk_size,
                            fmean(ratios),
                            incremental_mean,
                            incremental_mean / (duration / chunk_size),
                        )
                    )
    atomic_write(output_directory / SUMMARY_FILENAME, buffer.getvalue())


def materialize_results(
    plan: Sequence[Job],
    records: dict[str, sqlite3.Row],
    output_directory: Path,
) -> None:
    materialize_rows(plan, records, output_directory)
    materialize_summary(plan, records, output_directory)


def configure_and_build(
    stages: Sequence[str],
    build_directory: Path,
    build_jobs: int,
    no_build: bool,
) -> None:
    targets = [f"sttt_online_{stage.replace('/', '_')}" for stage in stages]
    if not targets:
        return
    binaries = [build_directory / "bin" / target for target in targets]
    if no_build:
        missing = [str(path) for path in binaries if not path.is_file()]
        if missing:
            raise RuntimeError(
                "missing online binaries: "
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
            *targets,
        ],
        cwd=PROJECT_ROOT,
        check=True,
    )


def print_targets(stages: Sequence[str]) -> None:
    for stage in stages:
        print(f"sttt_online_{stage.replace('/', '_')}")


def print_plan(
    args: argparse.Namespace,
    plan: Sequence[Job],
    stages: Sequence[str],
) -> None:
    selected = [job for job in plan if job.stage in stages]
    counts = Counter(job.stage for job in selected)
    print("Online experiment plan:")
    for stage in stages:
        print(f"  {stage}: {counts[stage]} configurations")
    print(f"Total: {len(selected)} configurations")
    print(f"Repetitions per configuration: {args.repetitions}")
    print("Online C++ monitors have no timeout.")


def initialize_online_checkpoint(
    connection: sqlite3.Connection,
    fingerprint: str,
    source: str,
    config: str,
) -> bool:
    """Initialize metadata or adopt changed inputs while preserving job rows."""
    metadata = dict(connection.execute("SELECT key, value FROM metadata"))
    previous = metadata.get("fingerprint")
    if previous is None:
        initialize_checkpoint(connection, fingerprint, source, config)
        return False
    if previous == fingerprint:
        return False

    with connection:
        connection.executemany(
            """
            INSERT INTO metadata(key, value) VALUES (?, ?)
            ON CONFLICT(key) DO UPDATE SET value=excluded.value
            """,
            (
                ("fingerprint", fingerprint),
                ("source_digest", source),
                ("config", config),
            ),
        )
    return True


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
        connection = sqlite3.connect(f"file:{state_path}?mode=ro", uri=True)
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

    selected = [job for job in plan if job.stage in stages]
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
            f"{counts['error']} failed, {not_started} not started"
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
        f"{counts['error']} failed, {not_started} not started."
    )
    print(f"Results: {output_directory}")


def _run(args: argparse.Namespace) -> int:
    stages = selected_stages(args)
    plan = build_plan(args)

    if args.list:
        print_targets(stages)
        return 0
    if args.dry_run:
        print_plan(args, plan, stages)
        return 0

    state_path = args.output_dir / STATE_FILENAME
    if args.status:
        return status_report(state_path, plan, stages, args)

    selected = [job for job in plan if job.stage in stages]
    paths = output_paths(plan, args.output_dir)
    if args.restart:
        reset_outputs(state_path, paths)
    else:
        ensure_fresh_outputs(state_path, paths)

    fingerprint, source, config = experiment_fingerprint(args)
    connection = connect_checkpoint(state_path)
    try:
        inputs_changed = initialize_online_checkpoint(
            connection, fingerprint, source, config
        )
        if inputs_changed:
            print(
                "Warning: checkpoint sources or settings changed; "
                "preserving completed configurations and continuing."
            )
        records = load_records(connection)
        materialize_results(plan, records, args.output_dir)

        pending = [
            job
            for job in selected
            if job.job_id not in records
            or records[job.job_id]["status"] != "complete"
        ]
        pending_stages = tuple(
            stage
            for stage in stages
            if any(job.stage == stage for job in pending)
        )
        configure_and_build(
            pending_stages,
            args.build_dir,
            args.build_jobs,
            args.no_build,
        )

        finished = sum(
            job.job_id in records and records[job.job_id]["status"] == "complete"
            for job in selected
        )
        print(
            f"Running {len(selected)} checkpointed configurations; "
            f"{finished} already finished."
        )

        previous_stage: str | None = None
        for job in selected:
            record = records.get(job.job_id)
            if record is not None and record["status"] == "complete":
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
                print("\nInterrupted; completed configurations are checkpointed.")
                print_summary(selected, records, args.output_dir)
                return 130

            records[job.job_id] = record_attempt(connection, job, attempt)
            if attempt.status == "complete":
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
                    (4,),
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
                print("The failed configuration will retry on the next invocation.")
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

    lock = acquire_lock(args.output_dir, LOCK_FILENAME, "online")
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
