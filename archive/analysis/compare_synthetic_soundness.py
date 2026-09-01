#!/usr/bin/env python3
"""Compare the active approximate and exact synthetic monitors off-suite.

The comparison checks the soundness inclusion, not verdict equality.  If the
exact monitor finds a violating trace, the approximate language must contain a
false behavior; if exact finds a satisfying trace, approximate must contain a
true behavior.  Consequently an exact-inconclusive result must also be
approximate-inconclusive, while an exact-conclusive result may be weakened to
inconclusive by the abstraction.
"""

from __future__ import annotations

import argparse
import importlib.util
import itertools
import subprocess
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType


PROJECT_ROOT = Path(__file__).resolve().parents[2]
EXACT_RUNNER_PATH = PROJECT_ROOT / "benchmarks/offline/exact/run.py"
APPROXIMATE_SOURCE = Path(__file__).with_name(
    "approximate_synthetic_query.cpp"
)
DEFAULT_HELPER = PROJECT_ROOT / "build/analysis/approximate_synthetic_query"
FORMULAS = ("phi1", "phi2", "phi3", "phi4", "phi5", "phi6")


@dataclass(frozen=True)
class Query:
    identifier: int
    family: str
    formula: str
    epsilon: int
    p_word: str
    q_word: str
    description: str


@dataclass(frozen=True)
class ApproximateResult:
    segments: int
    false_possible: bool
    true_possible: bool
    verdict: int


@dataclass(frozen=True)
class ExactResult:
    proves_true: bool
    proves_false: bool
    verdict: int

    @property
    def false_exists(self) -> bool:
        return not self.proves_true

    @property
    def true_exists(self) -> bool:
        return not self.proves_false


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mode",
        choices=("exhaustive", "shipped", "both"),
        default="exhaustive",
    )
    parser.add_argument(
        "--formula",
        action="append",
        choices=FORMULAS,
        dest="formulas",
        help="formula to compare; repeat for several (default: all)",
    )
    parser.add_argument(
        "--length",
        action="append",
        type=int,
        dest="lengths",
        help="exhaustive Boolean-word length (default: 4)",
    )
    parser.add_argument(
        "--duration",
        action="append",
        type=int,
        dest="durations",
        help="shipped trace duration (default: 4)",
    )
    parser.add_argument(
        "--epsilon",
        action="append",
        type=int,
        dest="epsilons",
        help="clock-skew value; repeat for several (default: 1,2)",
    )
    parser.add_argument(
        "--sample-limit",
        type=int,
        default=10,
        help="number of shipped samples starting at zero (default: 10)",
    )
    parser.add_argument(
        "--helper",
        type=Path,
        default=DEFAULT_HELPER,
        help="path for the standalone approximate helper",
    )
    parser.add_argument(
        "--no-build",
        action="store_true",
        help="reuse --helper without compiling it",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=250,
        help="print exact-solver progress every N comparisons; zero disables",
    )
    parser.add_argument(
        "--max-reports",
        type=int,
        default=20,
        help="maximum detailed violations/errors to print",
    )
    args = parser.parse_args()

    args.formulas = tuple(args.formulas or FORMULAS)
    args.lengths = tuple(args.lengths or (4,))
    args.durations = tuple(args.durations or (4,))
    args.epsilons = tuple(args.epsilons or (1, 2))

    if any(length < 1 for length in args.lengths):
        parser.error("--length must be positive")
    if any(duration < 1 for duration in args.durations):
        parser.error("--duration must be positive")
    if any(epsilon < 1 for epsilon in args.epsilons):
        parser.error("--epsilon must be positive")
    if not 1 <= args.sample_limit <= 100:
        parser.error("--sample-limit must lie between 1 and 100")
    if args.progress_every < 0:
        parser.error("--progress-every must not be negative")
    if args.max_reports < 0:
        parser.error("--max-reports must not be negative")
    return args


def load_exact_runner() -> ModuleType:
    specification = importlib.util.spec_from_file_location(
        "sttt_soundness_exact_runner", EXACT_RUNNER_PATH
    )
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load exact runner from {EXACT_RUNNER_PATH}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def compile_helper(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    command = [
        "g++",
        "-std=c++17",
        "-O2",
        "-I",
        str(PROJECT_ROOT / "include"),
        str(APPROXIMATE_SOURCE),
        "-o",
        str(path),
    ]
    subprocess.run(command, cwd=PROJECT_ROOT, check=True)


def all_words(length: int) -> tuple[str, ...]:
    return tuple(
        "".join(bits) for bits in itertools.product("01", repeat=length)
    )


def exhaustive_queries(args: argparse.Namespace, start: int) -> list[Query]:
    queries: list[Query] = []
    identifier = start
    for length in args.lengths:
        words = all_words(length)
        for epsilon in args.epsilons:
            if epsilon > length:
                continue
            for formula in args.formulas:
                for p_word in words:
                    for q_word in words:
                        queries.append(
                            Query(
                                identifier,
                                "exhaustive",
                                formula,
                                epsilon,
                                p_word,
                                q_word,
                                f"length={length}",
                            )
                        )
                        identifier += 1
    return queries


def boolean_word(rows: list[list[float]], duration: int) -> str:
    if len(rows) < duration:
        raise ValueError("shipped signal is shorter than its duration")
    for index, row in enumerate(rows[:duration]):
        if row[0] != index:
            raise ValueError(
                "shipped synthetic timestamps must be the integer grid"
            )
    return "".join("1" if row[1] > 0 else "0" for row in rows[:duration])


def shipped_queries(
    args: argparse.Namespace, runner: ModuleType, start: int
) -> list[Query]:
    queries: list[Query] = []
    identifier = start
    for duration in args.durations:
        for epsilon in args.epsilons:
            if epsilon > duration:
                continue
            for sample in range(args.sample_limit):
                p_word = boolean_word(
                    runner.random_signal(duration, sample), duration
                )
                q_word = boolean_word(
                    runner.random_signal(duration, sample + 100), duration
                )
                for formula in args.formulas:
                    queries.append(
                        Query(
                            identifier,
                            "shipped",
                            formula,
                            epsilon,
                            p_word,
                            q_word,
                            f"duration={duration} sample={sample}",
                        )
                    )
                    identifier += 1
    return queries


def run_approximate(
    helper: Path, queries: list[Query]
) -> dict[int, ApproximateResult]:
    input_text = "".join(
        f"{query.identifier} {query.formula} {query.epsilon} "
        f"{query.p_word} {query.q_word}\n"
        for query in queries
    )
    completed = subprocess.run(
        [str(helper)],
        input=input_text,
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "approximate helper failed:\n" + completed.stderr.rstrip()
        )

    results: dict[int, ApproximateResult] = {}
    for line in completed.stdout.splitlines():
        fields = line.split()
        if len(fields) != 5:
            raise RuntimeError(f"malformed approximate output: {line}")
        identifier, segments, can_false, can_true, verdict = map(int, fields)
        results[identifier] = ApproximateResult(
            segments,
            bool(can_false),
            bool(can_true),
            verdict,
        )
    if len(results) != len(queries):
        raise RuntimeError(
            f"approximate helper returned {len(results)} of {len(queries)} rows"
        )
    return results


def trace(word: str) -> list[list[float]]:
    rows = [
        [float(timestamp), float(int(value))]
        for timestamp, value in enumerate(word)
    ]
    rows.append([float(len(word)), float(int(word[-1]))])
    return rows


def run_exact(
    runner: ModuleType,
    monitor: ModuleType,
    query: Query,
) -> ExactResult:
    p_data = trace(query.p_word)
    q_data = trace(query.q_word)
    proves_true = runner.evaluate_random_check(
        query.formula,
        False,
        query.epsilon,
        p_data,
        q_data,
        monitor,
    )
    proves_false = runner.evaluate_random_check(
        query.formula,
        True,
        query.epsilon,
        p_data,
        q_data,
        monitor,
    )
    return ExactResult(
        proves_true,
        proves_false,
        runner.exact_verdict(proves_true, proves_false),
    )


def soundness_errors(
    exact: ExactResult, approximate: ApproximateResult
) -> list[str]:
    errors: list[str] = []
    if exact.proves_true and exact.proves_false:
        errors.append("exact encoding proves both the formula and its negation")
    if not approximate.false_possible and not approximate.true_possible:
        errors.append("approximate language contains no Boolean outcome")
    if exact.false_exists and not approximate.false_possible:
        errors.append("approximate language omits an exact violating trace")
    if exact.true_exists and not approximate.true_possible:
        errors.append("approximate language omits an exact satisfying trace")
    return errors


def print_summary(
    queries: list[Query],
    pair_counts: dict[tuple[str, str], Counter[tuple[int, int]]],
    violations: list[str],
    failures: list[str],
    elapsed: float,
) -> None:
    print(f"comparisons: {len(queries)}")
    print(f"exact seconds: {elapsed:.3f}")
    print(f"soundness violations: {len(violations)}")
    print(f"solver/analysis failures: {len(failures)}")
    print("verdict pairs (exact->approximate):")
    for family, formula in sorted(pair_counts):
        counts = pair_counts[(family, formula)]
        rendered = " ".join(
            f"{exact}->{approximate}:{count}"
            for (exact, approximate), count in sorted(counts.items())
        )
        print(f"  {family:10s} {formula}: {rendered}")


def main() -> int:
    args = parse_args()
    runner = load_exact_runner()
    monitor = runner.load_module("random")

    if not args.no_build:
        compile_helper(args.helper)
    elif not args.helper.is_file():
        raise FileNotFoundError(args.helper)

    queries: list[Query] = []
    if args.mode in ("exhaustive", "both"):
        queries.extend(exhaustive_queries(args, len(queries)))
    if args.mode in ("shipped", "both"):
        queries.extend(shipped_queries(args, runner, len(queries)))
    if not queries:
        raise RuntimeError("the selected matrix contains no comparisons")

    print(f"prepared {len(queries)} comparisons", flush=True)
    approximate_results = run_approximate(args.helper, queries)

    pair_counts: dict[tuple[str, str], Counter[tuple[int, int]]] = defaultdict(
        Counter
    )
    violations: list[str] = []
    failures: list[str] = []
    started = time.perf_counter()
    for completed, query in enumerate(queries, start=1):
        approximate = approximate_results[query.identifier]
        try:
            exact = run_exact(runner, monitor, query)
            pair_counts[(query.family, query.formula)][
                (exact.verdict, approximate.verdict)
            ] += 1
            errors = soundness_errors(exact, approximate)
            if errors:
                violations.append(
                    f"{query.family} {query.formula} eps={query.epsilon} "
                    f"p={query.p_word} q={query.q_word} {query.description}: "
                    f"exact={exact.verdict} "
                    f"(proves_true={int(exact.proves_true)},"
                    f"proves_false={int(exact.proves_false)}) "
                    f"approx={approximate.verdict} "
                    f"(false_possible={int(approximate.false_possible)},"
                    f"true_possible={int(approximate.true_possible)}) "
                    f"segments={approximate.segments}: {'; '.join(errors)}"
                )
        except Exception as error:  # keep the matrix running for diagnostics
            failures.append(
                f"{query.family} {query.formula} eps={query.epsilon} "
                f"p={query.p_word} q={query.q_word} {query.description}: {error}"
            )

        if args.progress_every and completed % args.progress_every == 0:
            print(
                f"checked {completed}/{len(queries)} "
                f"violations={len(violations)} failures={len(failures)}",
                flush=True,
            )

    elapsed = time.perf_counter() - started
    print_summary(queries, pair_counts, violations, failures, elapsed)
    for heading, entries in (
        ("violations", violations),
        ("failures", failures),
    ):
        if entries:
            print(f"{heading} (first {args.max_reports}):")
            for entry in entries[: args.max_reports]:
                print(f"  {entry}")
    return 1 if violations or failures else 0


if __name__ == "__main__":
    sys.exit(main())
