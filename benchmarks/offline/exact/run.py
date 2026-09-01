#!/usr/bin/env python3
"""Run the exact SMT offline monitoring experiments."""

from __future__ import annotations

import argparse
import importlib.util
import sys
import time
from pathlib import Path
from types import ModuleType
from typing import Callable, NamedTuple, Sequence

PROJECT_ROOT = Path(__file__).resolve().parents[3]
SMT_DIRECTORY = Path(__file__).resolve().parent / "smt"
SMT_MODULE_PREFIX = "_sttt_exact_"
DEFAULT_OUTPUT_DIRECTORY = PROJECT_ROOT / "results/offline/exact"
CASE_STUDY_MULTIPLIER = 20
CASE_STUDY_DURATION = 20


class Formula(NamedTuple):
    description: str
    positive_monitor: str
    negative_monitor: str
    negate_negative_inputs: bool = False
    arguments: tuple[int, ...] = ()


FORMULAS = {
    "phi1": Formula(
        "always (p and q)",
        "prog_always_conjunction",
        "prog_eventually_disjunction",
        negate_negative_inputs=True,
    ),
    "phi2": Formula(
        "always (p or q)",
        "prog_always_disjunction",
        "prog_eventually_conjunction",
        negate_negative_inputs=True,
    ),
    "phi3": Formula("p until q", "prog_until", "prog_not_until"),
    "phi4": Formula(
        "always (p implies eventually q)",
        "prog_always_implies_eventually",
        "prog_not_always_implies_eventually",
    ),
    "phi5": Formula(
        "always (p implies eventually_[0,1) q)",
        "prog_always_implies_eventually_timed",
        "prog_not_always_implies_eventually_timed",
        arguments=(0, 1),
    ),
    "phi6": Formula(
        "always (p implies eventually_[0,2) q)",
        "prog_always_implies_eventually_timed",
        "prog_not_always_implies_eventually_timed",
        arguments=(0, 2),
    ),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "subject",
        nargs="?",
        choices=("random", "mutual-separation", "water-tanks"),
        default="random",
    )
    parser.add_argument(
        "--formula",
        choices=(*FORMULAS, "all"),
        default="phi1",
        help="random-trace formula (default: phi1)",
    )
    parser.add_argument("--duration", type=int, default=8, help="random trace duration")
    parser.add_argument("--epsilon", type=int, default=1, help="clock-skew value")
    parser.add_argument("--sample", type=int, default=0, help="random sample index")
    parser.add_argument("--agents", type=int, choices=(2, 3, 4), default=2)
    parser.add_argument("--repeat", type=int, default=1)
    parser.add_argument(
        "--timeout",
        type=int,
        help="timeout per solver check in seconds",
    )
    parser.add_argument(
        "--full",
        action="store_true",
        help="run the reported parameter matrix instead of one instance",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIRECTORY,
    )
    parser.add_argument("--quiet", action="store_true")
    parser.add_argument(
        "--list",
        action="store_true",
        help="list subjects and formulas without importing Z3",
    )
    args = parser.parse_args()

    if args.repeat < 1:
        parser.error("--repeat must be positive")
    if args.timeout is not None and args.timeout < 1:
        parser.error("--timeout must be positive")
    if not 0 <= args.sample <= 99:
        parser.error("--sample must be between 0 and 99")
    return args


def load_module(name: str) -> ModuleType:
    module_name = f"{SMT_MODULE_PREFIX}{name}"
    if module_name in sys.modules:
        return sys.modules[module_name]

    path = SMT_DIRECTORY / f"{name}.py"
    specification = importlib.util.spec_from_file_location(module_name, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load SMT module from {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[module_name] = module

    try:
        specification.loader.exec_module(module)
    except ModuleNotFoundError as error:
        sys.modules.pop(module_name, None)
        if error.name == "z3":
            raise SystemExit(
                "Z3 is not installed. Install requirements-smt.txt, then retry."
            ) from error
        raise
    except BaseException:
        sys.modules.pop(module_name, None)
        raise
    return module


def read_columns(path: Path, count: int) -> list[list[float]]:
    rows: list[list[float]] = []
    with path.open(encoding="utf-8") as source:
        for line_number, line in enumerate(source, start=1):
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) < count:
                raise ValueError(f"{path}:{line_number}: expected {count} columns")
            rows.append([float(field) for field in fields[:count]])
    return rows


def random_signal(duration: int, sample: int) -> list[list[float]]:
    path = PROJECT_ROOT / "data/signals" / f"{duration}_{sample}.txt"
    rows = read_columns(path, 2)
    if len(rows) < duration:
        raise ValueError(f"{path}: expected at least {duration} samples")
    return rows


def exact_verdict(positive: bool, negative: bool) -> int:
    if positive and not negative:
        return 1
    if negative and not positive:
        return 0
    return 2


def evaluate_random_check(
    formula: str,
    negative: bool,
    epsilon: int,
    data_0: list[list[float]],
    data_1: list[list[float]],
    monitor: ModuleType,
) -> bool:
    try:
        specification = FORMULAS[formula]
    except KeyError as error:
        raise ValueError(f"unsupported formula: {formula}") from error

    if negative and specification.negate_negative_inputs:
        data_0 = monitor.negate(data_0)
        data_1 = monitor.negate(data_1)

    monitor_name = (
        specification.negative_monitor
        if negative
        else specification.positive_monitor
    )
    check: Callable[..., bool] = getattr(monitor, monitor_name)
    return check(epsilon, 1, data_0, data_1, *specification.arguments)


def run_random(args: argparse.Namespace) -> None:
    formulas = list(FORMULAS) if args.formula == "all" else [args.formula]
    monitor = load_module("random")
    monitor.set_solver_timeout(getattr(args, "timeout", None))

    durations = (4, 8, 16, 32) if args.full else (args.duration,)
    epsilons = (1, 2, 4, 8) if args.full else (args.epsilon,)
    samples: Sequence[int] = range(100) if args.full else (args.sample,)
    output_directory = args.output_dir / "random"
    output_directory.mkdir(parents=True, exist_ok=True)

    for formula in formulas:
        output_path = output_directory / f"{formula}.txt"
        with output_path.open("w", encoding="utf-8") as results:
            for duration in durations:
                for epsilon in epsilons:
                    if epsilon > duration:
                        continue
                    for sample in samples:
                        source_0 = random_signal(duration, sample)
                        source_1 = random_signal(duration, sample + 100)

                        preparation_seconds = 0.0
                        positive_seconds = 0.0
                        negative_seconds = 0.0
                        positive = False
                        negative = False

                        for _ in range(args.repeat):
                            started = time.perf_counter()
                            data_0 = monitor.preprocess(source_0, duration)
                            data_1 = monitor.preprocess(source_1, duration)
                            prepared = time.perf_counter()

                            positive = evaluate_random_check(
                                formula,
                                False,
                                epsilon,
                                data_0,
                                data_1,
                                monitor,
                            )
                            positive_finished = time.perf_counter()
                            negative = evaluate_random_check(
                                formula,
                                True,
                                epsilon,
                                data_0,
                                data_1,
                                monitor,
                            )
                            finished = time.perf_counter()

                            preparation_seconds += prepared - started
                            positive_seconds += positive_finished - prepared
                            negative_seconds += finished - positive_finished

                        prep = preparation_seconds / args.repeat
                        positive_time = prep + positive_seconds / args.repeat
                        negative_time = prep + negative_seconds / args.repeat
                        row = (
                            f"{duration} {epsilon} - {sample} - "
                            f"{positive_time} {negative_time} "
                            f"{int(positive)} {int(negative)} "
                            f"{exact_verdict(positive, negative)}"
                        )
                        results.write(row + "\n")
                        if not args.quiet:
                            print(formula, row, flush=True)


def scaled_case_data(
    subject: str, agents: int, multiplier: int, duration: int
) -> list[list[list[float]]]:
    if subject == "mutual-separation":
        directory = PROJECT_ROOT / "data/case_studies/uav"
        prefix = "s1_uav_"
        columns = 4
    else:
        directory = PROJECT_ROOT / "data/case_studies/water_tanks"
        prefix = "s5_tank_"
        columns = 2

    traces = []
    for agent in range(agents):
        path = directory / f"{prefix}{agent}"
        rows = read_columns(path, columns)
        if not rows:
            raise ValueError(f"{path}: expected at least one sample")
        for row in rows:
            scaled_timestamp = row[0] * multiplier
            rounded_timestamp = round(scaled_timestamp)
            if abs(scaled_timestamp - rounded_timestamp) > 1e-9:
                raise ValueError(
                    f"case-study timestamp {row[0]} does not map to the integer grid"
                )
            row[0] = rounded_timestamp
        if rows[-1][0] > duration:
            raise ValueError("case-study data extends beyond its declared duration")
        if rows[-1][0] < duration:
            rows.append([duration, *rows[-1][1:]])
        traces.append(rows)
    return traces


def run_case_study(args: argparse.Namespace) -> None:
    module = load_module("case_studies")
    module.set_solver_timeout(getattr(args, "timeout", None))
    if args.subject == "mutual-separation":
        positive_monitor = module.prove_mutual_separation
        negative_monitor = module.prove_mutual_separation_negative
        default_epsilons = (1, 2, 3, 4, 5)
        output_name = "mutual_separation.txt"
    else:
        positive_monitor = module.prove_water_tank
        negative_monitor = module.prove_water_tank_negative
        default_epsilons = (1, 2, 4, 8)
        output_name = "water_tanks.txt"

    agents_values = (2, 3, 4) if args.full else (args.agents,)
    epsilons = default_epsilons if args.full else (args.epsilon,)
    output_directory = args.output_dir / "case_studies"
    output_directory.mkdir(parents=True, exist_ok=True)

    with (output_directory / output_name).open("w", encoding="utf-8") as results:
        for agents in agents_values:
            traces = scaled_case_data(
                args.subject,
                agents,
                CASE_STUDY_MULTIPLIER,
                CASE_STUDY_DURATION,
            )
            for epsilon in epsilons:
                positive_elapsed = 0.0
                negative_elapsed = 0.0
                positive = False
                negative = False
                for _ in range(args.repeat):
                    started = time.perf_counter()
                    positive = positive_monitor(epsilon, traces)
                    positive_finished = time.perf_counter()
                    negative = negative_monitor(epsilon, traces)
                    finished = time.perf_counter()
                    positive_elapsed += positive_finished - started
                    negative_elapsed += finished - positive_finished

                row = (
                    f"1 {agents} {epsilon} "
                    f"{positive_elapsed / args.repeat} {int(positive)} "
                    f"{negative_elapsed / args.repeat} {int(negative)} "
                    f"{exact_verdict(positive, negative)}"
                )
                results.write(row + "\n")
                if not args.quiet:
                    print(args.subject, row, flush=True)


def main() -> None:
    args = parse_args()
    if args.list:
        print("subjects: random, mutual-separation, water-tanks")
        for formula, specification in FORMULAS.items():
            print(f"{formula}: {specification.description}")
        return

    try:
        if args.subject == "random":
            run_random(args)
        else:
            run_case_study(args)
    except TimeoutError as error:
        print(error, file=sys.stderr)
        raise SystemExit(124) from error


if __name__ == "__main__":
    main()
