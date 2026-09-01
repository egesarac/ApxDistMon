#!/usr/bin/env python3
"""Regression tests for the exact offline case-study SMT monitors."""

from __future__ import annotations

import importlib.util
import unittest
from itertools import product
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest import mock

try:
    import z3
except ModuleNotFoundError:
    z3 = None


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RUNNER_PATH = PROJECT_ROOT / "benchmarks/offline/exact/run.py"
REGRESSION_ORACLE_PATH = (
    PROJECT_ROOT
    / "tests/oracles/property_specialized_cut_inverse_clock.py"
)


def load_path(name, path):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load exact SMT module from {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def load_runner():
    return load_path("sttt_exact_case_study_runner", RUNNER_PATH)


@unittest.skipUnless(z3 is not None, "z3-solver is not installed")
class ExactCaseStudyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.runner = load_runner()
        cls.case_studies = cls.runner.load_module("case_studies")
        cls.random_monitor = cls.runner.load_module("random")
        cls.specialized_oracle = load_path(
            "sttt_case_study_oracle",
            REGRESSION_ORACLE_PATH,
        )

    def prove_water(self, epsilon, *traces):
        return self.case_studies.prove_water_tank(epsilon, list(traces))

    def prove_water_negative(self, epsilon, *traces):
        return self.case_studies.prove_water_tank_negative(
            epsilon, list(traces)
        )

    def prove_separation(self, epsilon, *traces):
        return self.case_studies.prove_mutual_separation(
            epsilon, list(traces)
        )

    def prove_separation_negative(self, epsilon, *traces):
        return self.case_studies.prove_mutual_separation_negative(
            epsilon, list(traces)
        )

    def test_solver_timeout_is_applied_in_milliseconds(self) -> None:
        solver = mock.Mock()
        self.addCleanup(self.case_studies.set_solver_timeout, None)
        self.case_studies.set_solver_timeout(360)

        with mock.patch.object(
            self.case_studies, "Solver", return_value=solver
        ):
            self.assertIs(self.case_studies._new_solver(), solver)

        solver.set.assert_called_once_with(timeout=360_000)

    @staticmethod
    def scalar_trace(*values: float) -> list[list[float]]:
        rows = [[timestamp, value] for timestamp, value in enumerate(values)]
        rows.append([len(values), values[-1]])
        return rows

    @staticmethod
    def point_trace(*values: tuple[float, float, float]) -> list[list[float]]:
        rows = [
            [timestamp, *value] for timestamp, value in enumerate(values)
        ]
        rows.append([len(values), *values[-1]])
        return rows

    def test_water_tank_uses_closed_violation_boundary(self) -> None:
        cases = (
            (self.prove_water, self.prove_water_negative, (5, 5)),
            (self.prove_water, self.prove_water_negative, (2, 3, 5)),
            (self.prove_water, self.prove_water_negative, (1, 2, 3, 4)),
        )
        for positive_monitor, negative_monitor, values in cases:
            with self.subTest(agents=len(values)):
                traces = [self.scalar_trace(value) for value in values]
                self.assertFalse(positive_monitor(1, *traces))
                self.assertTrue(negative_monitor(1, *traces))

    def test_water_tank_safe_constants_are_proved(self) -> None:
        cases = (
            (self.prove_water, (6, 5)),
            (self.prove_water, (4, 4, 4)),
            (self.prove_water, (3, 3, 3, 3)),
        )
        for monitor, values in cases:
            with self.subTest(agents=len(values)):
                traces = [self.scalar_trace(value) for value in values]
                self.assertTrue(monitor(1, *traces))

    def test_both_water_polarities_construct_the_full_flow(self) -> None:
        traces = [self.scalar_trace(5), self.scalar_trace(5)]
        with mock.patch.object(
            self.case_studies,
            "_build_full_flow",
            wraps=self.case_studies._build_full_flow,
        ) as build_full_flow:
            self.assertTrue(self.prove_water_negative(1, *traces))
        build_full_flow.assert_called_once()

    def test_last_water_sample_is_not_treated_as_terminal(self) -> None:
        falling = self.scalar_trace(6, 0)
        self.assertFalse(self.prove_water(1, falling, falling))

    def test_water_tank_three_valued_verdicts(self) -> None:
        cases = (
            (self.scalar_trace(6), self.scalar_trace(5), (True, False, 1)),
            (self.scalar_trace(5), self.scalar_trace(5), (False, True, 0)),
            (
                self.scalar_trace(6, 0),
                self.scalar_trace(6, 12),
                (False, False, 2),
            ),
        )
        for first, second, expected in cases:
            with self.subTest(expected=expected):
                positive = self.prove_water(1, first, second)
                negative = self.prove_water_negative(1, first, second)
                self.assertEqual(
                    (
                        positive,
                        negative,
                        self.runner.exact_verdict(positive, negative),
                    ),
                    expected,
                )

    def test_boolean_water_tank_matches_random_exact_invariants(self) -> None:
        # For Boolean signals, sum > 1 is p and q, while sum > 0 is p or q.
        # This directly compares both proof polarities of the independently
        # implemented case-study and random exact encodings.
        formula_thresholds = (("phi1", 1), ("phi2", 0))
        for duration in range(1, 4):
            words = tuple(product((0, 1), repeat=duration))
            for first_values in words:
                for second_values in words:
                    first = self.scalar_trace(*first_values)
                    second = self.scalar_trace(*second_values)
                    for formula, threshold in formula_thresholds:
                        with self.subTest(
                            duration=duration,
                            first=first_values,
                            second=second_values,
                            formula=formula,
                        ):
                            random_flags = tuple(
                                self.runner.evaluate_random_check(
                                    formula,
                                    negative,
                                    1,
                                    first,
                                    second,
                                    self.random_monitor,
                                )
                                for negative in (False, True)
                            )
                            case_study_flags = (
                                self.case_studies.prove_water_tank(
                                    1, [first, second], threshold
                                ),
                                self.case_studies.prove_water_tank_negative(
                                    1, [first, second], threshold
                                ),
                            )
                            self.assertEqual(
                                case_study_flags, random_flags
                            )

    def test_boolean_separation_matches_random_exact_disjunction(self) -> None:
        # (p, 0, 0) and (0, q, 0) are distinct exactly when p or q.
        for duration in range(1, 4):
            words = tuple(product((0, 1), repeat=duration))
            for first_values in words:
                for second_values in words:
                    first = self.scalar_trace(*first_values)
                    second = self.scalar_trace(*second_values)
                    first_points = self.point_trace(
                        *((value, 0, 0) for value in first_values)
                    )
                    second_points = self.point_trace(
                        *((0, value, 0) for value in second_values)
                    )
                    with self.subTest(
                        duration=duration,
                        first=first_values,
                        second=second_values,
                    ):
                        random_flags = tuple(
                            self.runner.evaluate_random_check(
                                "phi2",
                                negative,
                                1,
                                first,
                                second,
                                self.random_monitor,
                            )
                            for negative in (False, True)
                        )
                        case_study_flags = (
                            self.prove_separation(
                                1, first_points, second_points
                            ),
                            self.prove_separation_negative(
                                1, first_points, second_points
                            ),
                        )
                        self.assertEqual(case_study_flags, random_flags)

    def test_full_flow_matches_archived_specialized_oracle(self) -> None:
        words = tuple(product((0, 1), repeat=2))
        for agent_count in (2, 3, 4):
            threshold = agent_count // 2
            for agent_words in product(words, repeat=agent_count):
                scalar_traces = [
                    self.scalar_trace(*word) for word in agent_words
                ]
                point_traces = [
                    self.point_trace(
                        *((agent + value, 0, 0) for value in word)
                    )
                    for agent, word in enumerate(agent_words)
                ]
                with self.subTest(
                    agents=agent_count,
                    words=agent_words,
                ):
                    full_flow_water = (
                        self.case_studies.prove_water_tank(
                            1, scalar_traces, threshold
                        ),
                        self.case_studies.prove_water_tank_negative(
                            1, scalar_traces, threshold
                        ),
                    )
                    specialized_water = (
                        self.specialized_oracle.prove_water_tank(
                            1, scalar_traces, threshold
                        ),
                        self.specialized_oracle.prove_water_tank_negative(
                            1, scalar_traces, threshold
                        ),
                    )
                    self.assertEqual(full_flow_water, specialized_water)

                    full_flow_separation = (
                        self.case_studies.prove_mutual_separation(
                            1, point_traces
                        ),
                        self.case_studies.prove_mutual_separation_negative(
                            1, point_traces
                        ),
                    )
                    specialized_separation = (
                        self.specialized_oracle.prove_mutual_separation(
                            1, point_traces
                        ),
                        self.specialized_oracle.prove_mutual_separation_negative(
                            1, point_traces
                        ),
                    )
                    self.assertEqual(
                        full_flow_separation,
                        specialized_separation,
                    )

    def test_variation_detection_uses_all_position_coordinates(self) -> None:
        data = (
            (0, 0, 0),
            (0, 1, 0),
            (0, 1, 1),
        )
        self.assertEqual(
            self.case_studies._variation_points(data),
            (0, 1, 2, 3),
        )

    def test_mutual_separation_checks_every_pair(self) -> None:
        origin = self.point_trace((0, 0, 0))
        far = self.point_trace((9, 9, 9))
        other = self.point_trace((4, 4, 4))

        self.assertFalse(self.prove_separation(1, origin, origin))
        self.assertTrue(self.prove_separation_negative(1, origin, origin))
        self.assertFalse(
            self.prove_separation(1, far, origin, origin)
        )
        self.assertTrue(
            self.prove_separation_negative(1, far, origin, origin)
        )
        self.assertFalse(
            self.prove_separation(1, far, other, origin, origin)
        )
        self.assertTrue(
            self.prove_separation_negative(1, far, other, origin, origin)
        )

    def test_mutual_separation_safe_constants_are_proved(self) -> None:
        points = (
            self.point_trace((0, 0, 0)),
            self.point_trace((1, 0, 0)),
            self.point_trace((2, 0, 0)),
            self.point_trace((3, 0, 0)),
        )
        self.assertTrue(self.prove_separation(1, *points[:2]))
        self.assertTrue(self.prove_separation(1, *points[:3]))
        self.assertTrue(self.prove_separation(1, *points))

    def test_last_position_sample_is_not_treated_as_terminal(self) -> None:
        first = self.point_trace((0, 0, 0), (1, 1, 1))
        second = self.point_trace((2, 2, 2), (1, 1, 1))
        self.assertFalse(self.prove_separation(1, first, second))

    def test_mutual_separation_three_valued_verdicts(self) -> None:
        cases = (
            (
                self.point_trace((0, 0, 0)),
                self.point_trace((1, 0, 0)),
                (True, False, 1),
            ),
            (
                self.point_trace((0, 0, 0)),
                self.point_trace((0, 0, 0)),
                (False, True, 0),
            ),
            (
                self.point_trace((0, 0, 0), (1, 0, 0)),
                self.point_trace((1, 0, 0), (2, 0, 0)),
                (False, False, 2),
            ),
        )
        for first, second, expected in cases:
            with self.subTest(expected=expected):
                positive = self.prove_separation(1, first, second)
                negative = self.prove_separation_negative(1, first, second)
                self.assertEqual(
                    (
                        positive,
                        negative,
                        self.runner.exact_verdict(positive, negative),
                    ),
                    expected,
                )

    def test_simultaneous_swap_can_be_the_only_safe_flow(self) -> None:
        first = self.point_trace((0, 0, 0), (1, 0, 0))
        second = self.point_trace((1, 0, 0), (0, 0, 0))
        self.assertFalse(self.prove_separation(1, first, second))
        self.assertFalse(self.prove_separation_negative(1, first, second))

    def test_strict_reference_window_and_epsilon_hb_boundary(self) -> None:
        allowed = self.case_studies._allowed_local_times
        self.assertEqual(allowed(0, 3, 2, 2), (0,))
        self.assertEqual(allowed(2, 3, 2, 2), (0, 1))
        self.assertEqual(allowed(4, 3, 2, 2), (1, 2))

        first = self.scalar_trace(0, 1, 1)
        second = self.scalar_trace(0, 0, 1)
        encoding = self.case_studies._build_full_flow(
            1, [first, second], 1
        )
        first_occurrence = self.case_studies._edge_occurrence(
            encoding.clocks[0],
            1,
            encoding.lower_tick,
            encoding.upper_tick,
        )
        second_occurrence = self.case_studies._edge_occurrence(
            encoding.clocks[1],
            2,
            encoding.lower_tick,
            encoding.upper_tick,
        )
        encoding.solver.add(first_occurrence >= second_occurrence)
        self.assertEqual(encoding.solver.check(), z3.unsat)

    def test_runner_provides_dense_terminal_marked_traces(self) -> None:
        water = self.runner.scaled_case_data("water-tanks", 2, 20, 20)
        separation = self.runner.scaled_case_data(
            "mutual-separation", 2, 20, 20
        )

        for traces in (water, separation):
            for trace in traces:
                self.assertEqual(
                    [row[0] for row in trace], list(range(21))
                )

    def test_runner_writes_both_flags_and_three_valued_verdict(self) -> None:
        with TemporaryDirectory() as output_directory:
            args = SimpleNamespace(
                subject="water-tanks",
                full=False,
                agents=2,
                epsilon=1,
                output_dir=Path(output_directory),
                repeat=1,
                quiet=True,
            )
            self.runner.run_case_study(args)
            output_path = (
                Path(output_directory) / "case_studies/water_tanks.txt"
            )
            fields = output_path.read_text(encoding="utf-8").split()

        self.assertEqual(len(fields), 8)
        self.assertEqual(fields[:3], ["1", "2", "1"])
        self.assertEqual((fields[4], fields[6], fields[7]), ("1", "0", "1"))

    def test_decimal_threshold_counterexamples(self) -> None:
        boundary_values = (3.3333336, 3.3333336, 3.3333328)
        boundary = [self.scalar_trace(value) for value in boundary_values]
        self.assertFalse(self.prove_water(1, *boundary))
        self.assertTrue(self.prove_water_negative(1, *boundary))

        barely_safe = [
            self.scalar_trace(5.0000004),
            self.scalar_trace(5.0000004),
        ]
        self.assertTrue(self.prove_water(1, *barely_safe))
        self.assertFalse(self.prove_water_negative(1, *barely_safe))

        origin = self.point_trace((0, 0, 0))
        nearby = self.point_trace((0.0001, 0, 0))
        self.assertTrue(self.prove_separation(1, origin, nearby))
        self.assertFalse(self.prove_separation_negative(1, origin, nearby))

    def test_unknown_is_not_reported_as_a_proof(self) -> None:
        scalar = self.scalar_trace(6)
        low_scalar = self.scalar_trace(0)
        point = self.point_trace((0, 0, 0))
        cases = (
            (self.prove_water, (scalar, scalar)),
            (self.prove_water, (scalar, scalar, scalar)),
            (self.prove_water, (scalar, scalar, scalar, scalar)),
            (self.prove_water_negative, (low_scalar, low_scalar)),
            (
                self.prove_water_negative,
                (low_scalar, low_scalar, low_scalar),
            ),
            (
                self.prove_water_negative,
                (low_scalar, low_scalar, low_scalar, low_scalar),
            ),
            (self.prove_separation, (point, point)),
            (self.prove_separation, (point, point, point)),
            (self.prove_separation, (point, point, point, point)),
            (self.prove_separation_negative, (point, point)),
            (self.prove_separation_negative, (point, point, point)),
            (self.prove_separation_negative, (point, point, point, point)),
        )

        with mock.patch.object(z3.Solver, "check", return_value=z3.unknown):
            for monitor, traces in cases:
                with self.subTest(function=monitor.__name__):
                    with self.assertRaisesRegex(RuntimeError, "(?i)unknown"):
                        monitor(1, *traces)

    def test_timeout_unknown_is_reported_as_a_timeout(self) -> None:
        solver = mock.Mock()
        solver.check.return_value = z3.unknown
        solver.reason_unknown.return_value = "timeout"

        with self.assertRaisesRegex(TimeoutError, "(?i)timed out"):
            self.case_studies._check_unsat(solver)

    def test_epsilon_must_be_a_positive_integer(self) -> None:
        scalar = self.scalar_trace(6)
        point = self.point_trace((0, 0, 0))
        cases = (
            (self.prove_water, (scalar, scalar)),
            (self.prove_water_negative, (scalar, scalar)),
            (self.prove_separation, (point, point)),
            (self.prove_separation_negative, (point, point)),
        )
        for epsilon in (0, -1, 1.0, True):
            for monitor, traces in cases:
                with self.subTest(epsilon=epsilon, function=monitor.__name__):
                    with self.assertRaises(ValueError):
                        monitor(epsilon, *traces)


if __name__ == "__main__":
    unittest.main()
