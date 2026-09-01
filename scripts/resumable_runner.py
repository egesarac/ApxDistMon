"""Shared checkpoint, locking, and result helpers for experiment runners."""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import os
import shlex
import shutil
import sqlite3
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence


def positive_integer(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be positive")
    return parsed


def format_attempt_timing(
    status: str,
    result_row: str | None,
    elapsed_seconds: float,
    timing_columns: Sequence[int],
    unit: str,
) -> str:
    """Format the exact maximum result timing or an incomplete wall time."""
    if status == "complete" and result_row is not None:
        fields = result_row.split()
        try:
            timings = tuple(fields[index] for index in timing_columns)
            maximum = max(timings, key=float)
        except (IndexError, ValueError) as error:
            raise RuntimeError(
                "completed result contains an invalid timing"
            ) from error
        return f"avg {maximum}{unit}"
    return f"wall {elapsed_seconds:.3f}s"


def _cmake_cache_values(build_directory: Path) -> dict[str, str]:
    cache_path = build_directory / "CMakeCache.txt"
    if not cache_path.is_file():
        return {}

    values: dict[str, str] = {}
    for line in cache_path.read_text(
        encoding="utf-8", errors="replace"
    ).splitlines():
        if not line or line.startswith(("//", "#")) or "=" not in line:
            continue
        key_and_type, value = line.split("=", 1)
        key = key_and_type.split(":", 1)[0]
        values[key] = value
    return values


def _resolved_compiler(cache: dict[str, str]) -> str:
    compiler = cache.get("CMAKE_CXX_COMPILER", "")
    if not compiler or compiler.endswith("-NOTFOUND"):
        configured = shlex.split(os.environ.get("CXX", "c++"))
        compiler = configured[0] if configured else "c++"

    located = shutil.which(compiler)
    if located is not None:
        return str(Path(located).resolve())

    candidate = Path(compiler).expanduser()
    if candidate.is_absolute():
        return str(candidate.resolve())
    return compiler


def _compiler_version(compiler: str) -> str:
    try:
        process = subprocess.run(
            [compiler, "--version"],
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError:
        return "unavailable"

    output = process.stdout or process.stderr
    for line in output.splitlines():
        if line.strip():
            return line.strip()
    return f"unavailable (exit {process.returncode})"


def build_context(build_directory: Path) -> dict[str, str]:
    """Return stable build metadata used to protect timing checkpoints."""
    directory = build_directory.expanduser().resolve()
    cache = _cmake_cache_values(directory)
    compiler = _resolved_compiler(cache)
    return {
        "build_directory": str(directory),
        "build_profile": "Release/O3",
        "cxx_compiler": compiler,
        "cxx_compiler_version": _compiler_version(compiler),
        "cxxflags": cache.get(
            "CMAKE_CXX_FLAGS", os.environ.get("CXXFLAGS", "")
        ),
    }


def digest_files(project_root: Path, paths: Iterable[Path]) -> str:
    digest = hashlib.sha256()
    ordered = sorted(
        {path.resolve() for path in paths},
        key=lambda path: str(path.relative_to(project_root)),
    )
    for path in ordered:
        relative = path.relative_to(project_root)
        digest.update(str(relative).encode("utf-8"))
        digest.update(b"\0")
        with path.open("rb") as source:
            while chunk := source.read(1024 * 1024):
                digest.update(chunk)
        digest.update(b"\0")
    return digest.hexdigest()


def experiment_fingerprint(
    project_root: Path,
    paths: Iterable[Path],
    config: object,
) -> tuple[str, str, str]:
    source = digest_files(project_root, paths)
    serialized = json.dumps(config, sort_keys=True, separators=(",", ":"))
    digest = hashlib.sha256(
        f"{source}\0{serialized}".encode("utf-8")
    ).hexdigest()
    return digest, source, serialized


def remove_checkpoint(path: Path) -> None:
    for candidate in (
        path,
        Path(f"{path}-journal"),
        Path(f"{path}-wal"),
        Path(f"{path}-shm"),
    ):
        candidate.unlink(missing_ok=True)


def reset_outputs(state_path: Path, paths: Iterable[Path]) -> None:
    remove_checkpoint(state_path)
    for path in set(paths):
        path.unlink(missing_ok=True)


def ensure_fresh_outputs(state_path: Path, paths: Iterable[Path]) -> None:
    if state_path.exists():
        return
    if any(path.exists() for path in paths):
        raise RuntimeError(
            "generated results exist without a matching checkpoint; "
            "use --restart or choose another --output-dir"
        )


def acquire_lock(
    output_directory: Path,
    lock_filename: str,
    runner_description: str,
):
    output_directory.mkdir(parents=True, exist_ok=True)
    lock = (output_directory / lock_filename).open("a+", encoding="utf-8")
    try:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError as error:
        lock.close()
        raise RuntimeError(
            f"another {runner_description} experiment runner is already active"
        ) from error
    lock.seek(0)
    lock.truncate()
    lock.write(f"{os.getpid()}\n")
    lock.flush()
    return lock


def connect_checkpoint(path: Path) -> sqlite3.Connection:
    path.parent.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(path)
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA journal_mode=WAL")
    connection.execute("PRAGMA synchronous=FULL")
    connection.executescript(
        """
        CREATE TABLE IF NOT EXISTS metadata (
            key TEXT PRIMARY KEY,
            value TEXT NOT NULL
        );
        CREATE TABLE IF NOT EXISTS jobs (
            job_id TEXT PRIMARY KEY,
            stage TEXT NOT NULL,
            status TEXT NOT NULL,
            result_row TEXT,
            elapsed_seconds REAL NOT NULL,
            message TEXT NOT NULL,
            finished_at TEXT NOT NULL
        );
        """
    )
    return connection


def initialize_checkpoint(
    connection: sqlite3.Connection,
    fingerprint: str,
    source: str,
    config: str,
) -> None:
    metadata = dict(connection.execute("SELECT key, value FROM metadata"))
    previous = metadata.get("fingerprint")
    if previous is not None and previous != fingerprint:
        raise RuntimeError(
            "the checkpoint belongs to different sources or settings; "
            "use --restart to discard it"
        )
    if previous is None:
        with connection:
            connection.executemany(
                "INSERT INTO metadata(key, value) VALUES (?, ?)",
                (
                    ("fingerprint", fingerprint),
                    ("source_digest", source),
                    ("config", config),
                    ("created_at", datetime.now(timezone.utc).isoformat()),
                ),
            )


def load_records(connection: sqlite3.Connection) -> dict[str, sqlite3.Row]:
    return {
        row["job_id"]: row
        for row in connection.execute(
            "SELECT job_id, stage, status, result_row, elapsed_seconds, "
            "message, finished_at FROM jobs"
        )
    }


def record_attempt(connection, job, attempt) -> sqlite3.Row:
    with connection:
        connection.execute(
            """
            INSERT INTO jobs(
                job_id, stage, status, result_row, elapsed_seconds,
                message, finished_at
            ) VALUES (?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(job_id) DO UPDATE SET
                stage=excluded.stage,
                status=excluded.status,
                result_row=excluded.result_row,
                elapsed_seconds=excluded.elapsed_seconds,
                message=excluded.message,
                finished_at=excluded.finished_at
            """,
            (
                job.job_id,
                job.stage,
                attempt.status,
                attempt.row,
                attempt.elapsed_seconds,
                attempt.message,
                datetime.now(timezone.utc).isoformat(),
            ),
        )
    record = connection.execute(
        "SELECT job_id, stage, status, result_row, elapsed_seconds, "
        "message, finished_at FROM jobs WHERE job_id = ?",
        (job.job_id,),
    ).fetchone()
    if record is None:
        raise RuntimeError(f"failed to checkpoint {job.job_id}")
    return record


def diagnostic_text(value: str | bytes | None) -> str:
    if value is None:
        return ""
    if isinstance(value, bytes):
        value = value.decode("utf-8", errors="replace")
    lines = [line.strip() for line in value.splitlines() if line.strip()]
    return " | ".join(lines[-4:])[-2000:]


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    with temporary.open("w", encoding="utf-8") as destination:
        destination.write(text)
        destination.flush()
        os.fsync(destination.fileno())
    os.replace(temporary, path)


def materialize_rows(
    plan: Sequence[object],
    records: dict[str, sqlite3.Row],
    output_directory: Path,
) -> None:
    materialized: dict[Path, list[str]] = {}
    for job in plan:
        rows = materialized.setdefault(job.output_path, [])
        record = records.get(job.job_id)
        if (
            record is not None
            and record["status"] == "complete"
            and record["result_row"] is not None
        ):
            rows.append(record["result_row"])

    for relative, rows in materialized.items():
        body = "".join(f"{row}\n" for row in rows)
        atomic_write(output_directory / relative, body)
