"""Exact continuous-clock checks for the offline safety case studies."""

from decimal import Decimal
from itertools import combinations, product
from math import isfinite
from numbers import Integral, Real as RealNumber

from z3 import And, If, Or, Real, RealVal, Solver, Sum, sat, unknown, unsat


_SOLVER_TIMEOUT_MILLISECONDS = None


def set_solver_timeout(seconds):
    global _SOLVER_TIMEOUT_MILLISECONDS
    if seconds is not None and (
        isinstance(seconds, bool)
        or not isinstance(seconds, Integral)
        or seconds <= 0
    ):
        raise ValueError("solver timeout must be a positive integer")
    _SOLVER_TIMEOUT_MILLISECONDS = (
        None if seconds is None else seconds * 1000
    )


def _new_solver():
    solver = Solver()
    if _SOLVER_TIMEOUT_MILLISECONDS is not None:
        solver.set(timeout=_SOLVER_TIMEOUT_MILLISECONDS)
    return solver


def _validate_epsilon(epsilon):
    if (
        isinstance(epsilon, bool)
        or not isinstance(epsilon, Integral)
        or epsilon <= 0
    ):
        raise ValueError("eps must be a positive integer")


def _validate_signals(signals, value_count):
    if not 2 <= len(signals) <= 4:
        raise ValueError("case studies require two, three, or four agents")

    normalized = []
    duration = None
    for agent, data in enumerate(signals):
        if len(data) < 2:
            raise ValueError(
                f"data_{agent} must contain a sample and a terminal marker"
            )

        rows = []
        previous_timestamp = None
        for index, row in enumerate(data):
            if len(row) < value_count + 1:
                raise ValueError(
                    f"data_{agent}[{index}] must contain a timestamp and "
                    f"{value_count} value(s)"
                )

            timestamp = row[0]
            if (
                isinstance(timestamp, bool)
                or not isinstance(timestamp, RealNumber)
                or not isfinite(float(timestamp))
                or int(timestamp) != timestamp
            ):
                raise ValueError("case-study timestamps must be integers")
            timestamp = int(timestamp)
            if previous_timestamp is not None and timestamp <= previous_timestamp:
                raise ValueError("case-study timestamps must be strictly increasing")
            if index == 0 and timestamp != 0:
                raise ValueError("case-study traces must start at timestamp zero")

            values = []
            for value in row[1:value_count + 1]:
                if (
                    isinstance(value, bool)
                    or not isinstance(value, RealNumber)
                    or not isfinite(float(value))
                ):
                    raise ValueError("case-study values must be finite real numbers")
                values.append(value)

            rows.append((timestamp, tuple(values)))
            previous_timestamp = timestamp

        signal_duration = rows[-1][0]
        if signal_duration <= 0:
            raise ValueError("case-study duration must be positive")
        if duration is None:
            duration = signal_duration
        elif duration != signal_duration:
            raise ValueError("all case-study traces must have the same duration")
        normalized.append(rows)

    return normalized, duration


def _piecewise_value(reading, rows, column):
    """Evaluate a right-continuous signal before its terminal marker."""
    value = RealVal(str(rows[0][1][column]))
    for timestamp, values in rows[1:-1]:
        value = If(reading >= timestamp, RealVal(str(values[column])), value)
    return value


def _state_model(epsilon, signals, value_count):
    """Encode one state of an admissible continuous consistent-cut flow."""
    _validate_epsilon(epsilon)
    epsilon = int(epsilon)
    rows, duration = _validate_signals(signals, value_count)
    solver = _new_solver()
    physical_time = Real("case_physical_time")
    readings = [
        Real(f"case_reading_{agent}") for agent in range(len(rows))
    ]

    # Every feasible knot extends through (0,0) and (d,d) by strict
    # piecewise-linear interpolation. Skew remains strict between the knots.
    solver.add(
        Or(
            And(physical_time == 0, *(reading == 0 for reading in readings)),
            And(
                0 < physical_time,
                physical_time < duration,
                *(And(0 < reading, reading < duration) for reading in readings),
            ),
        )
    )
    for reading in readings:
        solver.add(reading - physical_time < epsilon)
        solver.add(physical_time - reading < epsilon)

    for left, right in combinations(readings, 2):
        solver.add(left - right < epsilon, right - left < epsilon)

    values = [
        [
            _piecewise_value(readings[agent], data, column)
            for column in range(value_count)
        ]
        for agent, data in enumerate(rows)
    ]
    return solver, values


def _maximal_pieces(rows):
    """Return maximal half-open constant pieces, excluding the terminal row."""
    pieces = []
    for index, (start, values) in enumerate(rows[:-1]):
        end = rows[index + 1][0]
        if pieces and pieces[-1][2] == values:
            previous_start, _, previous_values = pieces[-1]
            pieces[-1] = (previous_start, end, previous_values)
        else:
            pieces.append((start, end, values))
    return pieces


def _safe_flow_model(epsilon, rows, duration):
    """Encode a complete admissible clock family by inverse-clock knots."""
    solver = _new_solver()
    agent_clocks = [
        [
            Real(f"case_inverse_{agent}_{local_time}")
            for local_time in range(duration + 1)
        ]
        for agent in range(len(rows))
    ]

    for clock in agent_clocks:
        solver.add(clock[0] == 0, clock[duration] == duration)
        solver.add(
            *(
                clock[local_time] < clock[local_time + 1]
                for local_time in range(duration)
            )
        )

    reference_clock = [
        RealVal(local_time) for local_time in range(duration + 1)
    ]
    clocks = [reference_clock, *agent_clocks]
    if epsilon <= duration:
        for target in range(len(clocks)):
            for source in range(len(clocks)):
                if target == source:
                    continue
                solver.add(
                    *(
                        clocks[source][local_time]
                        < clocks[target][local_time + epsilon]
                        for local_time in range(duration - epsilon + 1)
                    )
                )

    return solver, agent_clocks


def _forbid_common_intersection(solver, clocks, pieces):
    """Require the retimed half-open pieces to have empty intersection."""
    solver.add(
        Or(
            *(
                clocks[left][pieces[left][1]]
                <= clocks[right][pieces[right][0]]
                for left in range(len(pieces))
                for right in range(len(pieces))
                if left != right
            )
        )
    )


def _prove_no_counterexample(solver):
    result = solver.check()
    if result == unsat:
        return True
    if result == sat:
        return False
    if result == unknown:
        reason = solver.reason_unknown()
        if "timeout" in reason.lower():
            raise TimeoutError(f"Z3 timed out: {reason}")
        raise RuntimeError(f"Z3 returned unknown: {reason}")
    raise RuntimeError(f"unexpected Z3 result: {result}")


def prove_water_tank(epsilon, signals, threshold=10):
    """Prove G(sum(signals) > threshold) for every admissible flow."""
    solver, values = _state_model(epsilon, signals, 1)
    solver.add(Sum([value[0] for value in values]) <= RealVal(str(threshold)))
    return _prove_no_counterexample(solver)


def prove_mutual_separation(epsilon, signals):
    """Prove that every pair of three-dimensional positions stays distinct."""
    solver, values = _state_model(epsilon, signals, 3)
    collisions = [
        And(
            *(
                values[left][axis] == values[right][axis]
                for axis in range(3)
            )
        )
        for left, right in combinations(range(len(values)), 2)
    ]
    solver.add(Or(*collisions))
    return _prove_no_counterexample(solver)


def prove_water_tank_negative(epsilon, signals, threshold=10):
    """Prove that every admissible flow violates G(sum > threshold)."""
    _validate_epsilon(epsilon)
    epsilon = int(epsilon)
    rows, duration = _validate_signals(signals, 1)
    threshold_value = Decimal(str(threshold))
    piece_sets = [_maximal_pieces(data) for data in rows]

    # If even the pointwise minima exceed the threshold, no bad tuple exists.
    minimum_total = sum(
        min(Decimal(str(piece[2][0])) for piece in pieces)
        for pieces in piece_sets
    )
    if minimum_total > threshold_value:
        return False

    solver, clocks = _safe_flow_model(epsilon, rows, duration)
    for pieces in product(*piece_sets):
        total = sum(Decimal(str(piece[2][0])) for piece in pieces)
        if total <= threshold_value:
            _forbid_common_intersection(solver, clocks, pieces)

    return _prove_no_counterexample(solver)


def prove_mutual_separation_negative(epsilon, signals):
    """Prove that every admissible flow contains a drone collision."""
    _validate_epsilon(epsilon)
    epsilon = int(epsilon)
    rows, duration = _validate_signals(signals, 3)
    solver, clocks = _safe_flow_model(epsilon, rows, duration)
    piece_sets = [_maximal_pieces(data) for data in rows]

    for left, right in combinations(range(len(rows)), 2):
        for left_piece in piece_sets[left]:
            for right_piece in piece_sets[right]:
                if left_piece[2] == right_piece[2]:
                    _forbid_common_intersection(
                        solver,
                        [clocks[left], clocks[right]],
                        [left_piece, right_piece],
                    )

    return _prove_no_counterexample(solver)
