#!/usr/bin/env python3
"""Regression tests for the exact random-trace SMT monitors."""

from __future__ import annotations

import importlib.util
import random as stdlib_random
from functools import lru_cache
from itertools import combinations, product
import unittest
from pathlib import Path
from unittest import mock

try:
    import z3
except ModuleNotFoundError:
    z3 = None


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RUNNER_PATH = PROJECT_ROOT / "benchmarks/offline/exact/run.py"


def load_runner():
    specification = importlib.util.spec_from_file_location(
        "sttt_exact_random_runner", RUNNER_PATH
    )
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load exact SMT runner from {RUNNER_PATH}")
    runner = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(runner)
    return runner


@lru_cache(maxsize=None)
def _oracle_clock_paths(duration: int, epsilon: int):
    """Enumerate the finite clock paths independently of the SMT model."""
    scale = 2 * duration
    horizon = scale * duration
    paths = []

    for transition_ticks in combinations(
        range(1, horizon), duration - 1
    ):
        transitions = set(transition_ticks)
        reading = 0
        path = []
        for tick in range(horizon):
            if tick in transitions:
                reading += 1
            path.append(reading)

        admissible = True
        for tick, local_time in enumerate(path):
            if not any(
                local_tick // scale == local_time
                and abs(local_tick - tick) < scale * epsilon
                for local_tick in range(horizon)
            ):
                admissible = False
                break

        if admissible:
            paths.append(tuple(path))

    return tuple(paths)


def _oracle_variations(word):
    return tuple(
        index
        for index in range(1, len(word))
        if word[index] != word[index - 1]
    )


@lru_cache(maxsize=None)
def _oracle_admissible_pairs(word_0, word_1, epsilon):
    paths = _oracle_clock_paths(len(word_0), epsilon)
    variations_0 = _oracle_variations(word_0)
    variations_1 = _oracle_variations(word_1)
    variation_times = set(variations_0 + variations_1)
    occurrences = {
        path: {
            local_time: path.index(local_time)
            for local_time in variation_times
        }
        for path in paths
    }
    pairs = []

    for clock_0 in paths:
        for clock_1 in paths:
            if any(
                abs(reading_0 - reading_1) > epsilon
                for reading_0, reading_1 in zip(clock_0, clock_1)
            ):
                continue

            ordered = True
            for local_time_0 in variations_0:
                for local_time_1 in variations_1:
                    occurrence_0 = occurrences[clock_0][local_time_0]
                    occurrence_1 = occurrences[clock_1][local_time_1]
                    if (
                        local_time_0 + epsilon <= local_time_1
                        and occurrence_0 >= occurrence_1
                    ) or (
                        local_time_1 + epsilon <= local_time_0
                        and occurrence_1 >= occurrence_0
                    ):
                        ordered = False
                        break
                if not ordered:
                    break

            if ordered:
                pairs.append((clock_0, clock_1))

    if not pairs:
        raise AssertionError("the discrete clock model must be nonempty")
    return tuple(pairs)


def _oracle_response(clock_0, clock_1, word_0, word_1, epsilon, a, b):
    scale = 2 * len(word_0)
    horizon = len(clock_0)
    for tick in range(horizon):
        if not word_0[clock_0[tick]]:
            continue
        if not any(
            word_1[clock_1[witness]]
            for witness in range(
                tick + a * scale,
                min(tick + b * scale, horizon),
            )
        ):
            return False
    return True


def _oracle_verdict(word_0, word_1, epsilon, a, b):
    outcomes = tuple(
        _oracle_response(
            clock_0, clock_1, word_0, word_1, epsilon, a, b
        )
        for clock_0, clock_1 in _oracle_admissible_pairs(
            word_0, word_1, epsilon
        )
    )
    positive = all(outcomes)
    negative = all(not outcome for outcome in outcomes)
    verdict = 1 if positive else 0 if negative else 2
    return positive, negative, verdict


@unittest.skipUnless(z3 is not None, "z3-solver is not installed")
class RandomExactTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.runner = load_runner()
        cls.monitor = cls.runner.load_module("random")

    @staticmethod
    def trace(*values: int) -> list[list[float]]:
        """Build the piecewise-constant format produced by preprocess()."""
        samples = [
            [float(timestamp), float(value)]
            for timestamp, value in enumerate(values)
        ]
        samples.append([float(len(values)), float(values[-1])])
        return samples

    def evaluate(
        self,
        formula: str,
        data_0: list[list[float]],
        data_1: list[list[float]],
        epsilon: int = 1,
    ) -> tuple[bool, bool, int]:
        positive = self.runner.evaluate_random_check(
            formula, False, epsilon, data_0, data_1, self.monitor
        )
        negative = self.runner.evaluate_random_check(
            formula, True, epsilon, data_0, data_1, self.monitor
        )
        return positive, negative, self.runner.exact_verdict(positive, negative)

    def test_loader_distinguishes_exact_and_standard_library_random(self) -> None:
        exact_path = Path(self.monitor.__file__).resolve()
        stdlib_path = Path(stdlib_random.__file__).resolve()

        self.assertTrue(exact_path.as_posix().endswith("smt/random.py"))
        self.assertIsNot(self.monitor, stdlib_random)
        self.assertNotEqual(exact_path, stdlib_path)

    def test_solver_timeout_is_applied_in_milliseconds(self) -> None:
        solver = mock.Mock()
        self.addCleanup(self.monitor.set_solver_timeout, None)
        self.monitor.set_solver_timeout(30)

        with mock.patch.object(self.monitor, "Solver", return_value=solver):
            self.assertIs(self.monitor._new_solver(), solver)

        solver.set.assert_called_once_with(timeout=30_000)

    def test_constant_traces_cover_each_untimed_formula(self) -> None:
        cases = (
            ("phi1", self.trace(1, 1), self.trace(1, 1), (True, False, 1)),
            ("phi1", self.trace(0, 0), self.trace(1, 1), (False, True, 0)),
            ("phi2", self.trace(1, 1), self.trace(0, 0), (True, False, 1)),
            ("phi2", self.trace(0, 0), self.trace(0, 0), (False, True, 0)),
            ("phi3", self.trace(1, 1), self.trace(1, 1), (True, False, 1)),
            ("phi3", self.trace(1, 1), self.trace(0, 0), (False, True, 0)),
            ("phi4", self.trace(0, 0), self.trace(0, 0), (True, False, 1)),
            ("phi4", self.trace(1, 1), self.trace(0, 0), (False, True, 0)),
        )

        for formula, data_0, data_1, expected in cases:
            with self.subTest(formula=formula, expected=expected):
                self.assertEqual(self.evaluate(formula, data_0, data_1), expected)

    def test_until_excludes_the_q_witness_from_the_p_prefix(self) -> None:
        # q holds immediately, so p U q holds even though p is false there.
        self.assertEqual(
            self.evaluate("phi3", self.trace(0, 0), self.trace(1, 1)),
            (True, False, 1),
        )

    def test_until_keeps_both_admissible_edge_orders(self) -> None:
        # q-rise before p-fall satisfies until; p-fall before q-rise violates it.
        self.assertEqual(
            self.evaluate("phi3", self.trace(1, 0), self.trace(0, 1)),
            (False, False, 2),
        )

    def test_response_keeps_both_admissible_edge_orders(self) -> None:
        # p-fall first satisfies response; q-fall first leaves a pending p.
        self.assertEqual(
            self.evaluate("phi4", self.trace(1, 0), self.trace(1, 0)),
            (False, False, 2),
        )

    def test_unknown_is_not_reported_as_a_proof(self) -> None:
        data_0 = self.trace(1, 0)
        data_1 = self.trace(0, 1)
        monitors_and_arguments = (
            ("prog_not_always_implies_eventually", ()),
            ("prog_always_implies_eventually", ()),
            ("prog_always_conjunction", ()),
            ("prog_always_disjunction", ()),
            ("prog_eventually_conjunction", ()),
            ("prog_eventually_disjunction", ()),
            ("prog_until", ()),
            ("prog_not_until", ()),
        )

        with mock.patch.object(z3.Solver, "check", return_value=z3.unknown):
            for function_name, extra_arguments in monitors_and_arguments:
                with self.subTest(function=function_name):
                    monitor = getattr(self.monitor, function_name)
                    with self.assertRaisesRegex(RuntimeError, "(?i)unknown"):
                        monitor(1, 1, data_0, data_1, *extra_arguments)

    def test_timeout_unknown_is_reported_as_a_timeout(self) -> None:
        solver = mock.Mock()
        solver.check.return_value = z3.unknown
        solver.reason_unknown.return_value = "timeout"

        with self.assertRaisesRegex(TimeoutError, "(?i)timed out"):
            self.monitor._check_solver(solver)

    def test_epsilon_must_be_a_positive_integer(self) -> None:
        data_0 = self.trace(1, 0)
        data_1 = self.trace(0, 1)
        monitors_and_arguments = (
            ("prog_not_always_implies_eventually", ()),
            ("prog_always_implies_eventually", ()),
            ("prog_always_conjunction", ()),
            ("prog_always_disjunction", ()),
            ("prog_eventually_conjunction", ()),
            ("prog_eventually_disjunction", ()),
            ("prog_until", ()),
            ("prog_not_until", ()),
        )

        for epsilon in (0, -1, 1.0, True):
            for function_name, extra_arguments in monitors_and_arguments:
                with self.subTest(epsilon=epsilon, function=function_name):
                    monitor = getattr(self.monitor, function_name)
                    with self.assertRaises(ValueError):
                        monitor(
                            epsilon, 1, data_0, data_1, *extra_arguments
                        )


@unittest.skipUnless(z3 is not None, "z3-solver is not installed")
class TimedExactTests(unittest.TestCase):
    """Correctness checks for the standard timed exact baseline."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.runner = load_runner()
        cls.monitor = cls.runner.load_module("random")

    @staticmethod
    def trace(*values: int) -> list[list[float]]:
        samples = [
            [float(timestamp), float(value)]
            for timestamp, value in enumerate(values)
        ]
        samples.append([float(len(values)), float(values[-1])])
        return samples

    def evaluate_bounds(
        self,
        word_0: tuple[int, ...],
        word_1: tuple[int, ...],
        epsilon: int,
        a: int,
        b: int,
    ) -> tuple[bool, bool, int]:
        data_0 = self.trace(*word_0)
        data_1 = self.trace(*word_1)
        positive = (
            self.monitor.prog_always_implies_eventually_timed(
                epsilon, 1, data_0, data_1, a, b
            )
        )
        negative = (
            self.monitor.prog_not_always_implies_eventually_timed(
                epsilon, 1, data_0, data_1, a, b
            )
        )
        return (
            positive,
            negative,
            self.runner.exact_verdict(positive, negative),
        )

    def test_timed_baseline_matches_exhaustive_discrete_oracle(self) -> None:
        bounds = ((0, 1), (0, 2), (1, 2))
        for duration in range(1, 4):
            for word_0 in product((0, 1), repeat=duration):
                for word_1 in product((0, 1), repeat=duration):
                    for a, b in bounds:
                        with self.subTest(
                            duration=duration,
                            word_0=word_0,
                            word_1=word_1,
                            bounds=(a, b),
                        ):
                            actual = self.evaluate_bounds(
                                word_0, word_1, 1, a, b
                            )
                            expected = _oracle_verdict(
                                word_0, word_1, 1, a, b
                            )
                            self.assertEqual(actual, expected)
                            self.assertFalse(actual[0] and actual[1])

    def test_timed_baseline_matches_large_epsilon_oracle(self) -> None:
        cases = (
            ((0, 0, 1), (0, 1, 1), 0, 2),
            ((0, 1, 0), (0, 0, 1), 0, 1),
        )
        for word_0, word_1, a, b in cases:
            with self.subTest(word_0=word_0, word_1=word_1):
                self.assertEqual(
                    self.evaluate_bounds(word_0, word_1, 2, a, b),
                    _oracle_verdict(word_0, word_1, 2, a, b),
                )

    def test_timed_baseline_scales_grid_by_duration(self) -> None:
        encoding = mock.Mock()
        data_0 = self.trace(0, 1, 1, 0)
        data_1 = self.trace(1, 0, 0, 1)

        with (
            mock.patch.object(
                self.monitor,
                "_build_discrete_encoding",
                return_value=(encoding, False),
            ) as build,
            mock.patch.object(
                self.monitor, "_add_timed_response_query"
            ) as add_query,
            mock.patch.object(
                self.monitor, "_check_solver", return_value=z3.unsat
            ),
        ):
            self.assertTrue(
                self.monitor.prog_always_implies_eventually_timed(
                    1, 1, data_0, data_1, 0, 1
                )
            )

        self.assertEqual(build.call_args.kwargs["pad"], 8)
        add_query.assert_called_once_with(encoding, "response", 0, 1, 8)

    def test_timed_duration_scale_recovers_dense_counterexample(self) -> None:
        self.assertEqual(
            self.evaluate_bounds((0, 1, 1, 0), (1, 0, 0, 1), 1, 0, 1),
            (False, False, 2),
        )

    def test_timed_baseline_does_not_fabricate_a_tail_witness(self) -> None:
        self.assertEqual(
            self.evaluate_bounds((1,), (1,), 1, 1, 2),
            (False, True, 0),
        )

    def test_timed_baseline_rejects_invalid_parameters(self) -> None:
        data_0 = self.trace(1, 0)
        data_1 = self.trace(0, 1)
        monitor = self.monitor.prog_always_implies_eventually_timed

        for epsilon in (0, -1, 1.0, True):
            with self.subTest(epsilon=epsilon):
                with self.assertRaises(ValueError):
                    monitor(epsilon, 1, data_0, data_1, 0, 1)

        for bounds in ((0, 0), (-1, 1), (2, 1), (0.0, 1)):
            with self.subTest(bounds=bounds):
                with self.assertRaises(ValueError):
                    monitor(1, 1, data_0, data_1, *bounds)

    def test_timed_baseline_treats_segmentation_as_compatibility_only(self) -> None:
        data_0 = self.trace(0, 1, 0, 0)
        data_1 = self.trace(0, 0, 0, 1)
        for polarity in ("", "not_"):
            monitor = getattr(
                self.monitor,
                "prog_"
                + polarity
                + "always_implies_eventually_timed",
            )
            with self.subTest(polarity=polarity):
                self.assertEqual(
                    monitor(1, 1, data_0, data_1, 0, 1),
                    monitor(1, 4, data_0, data_1, 0, 1),
                )

    def test_timed_baseline_never_reports_unknown_as_a_proof(self) -> None:
        data_0 = self.trace(1, 0)
        data_1 = self.trace(0, 1)
        with mock.patch.object(z3.Solver, "check", return_value=z3.unknown):
            for polarity in ("", "not_"):
                monitor = getattr(
                    self.monitor,
                    "prog_"
                    + polarity
                    + "always_implies_eventually_timed",
                )
                with self.subTest(polarity=polarity):
                    with self.assertRaisesRegex(RuntimeError, "(?i)unknown"):
                        monitor(1, 1, data_0, data_1, 0, 1)


if __name__ == "__main__":
    unittest.main()
