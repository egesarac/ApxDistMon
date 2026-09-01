"""Corrected RV-derived full-flow SMT monitors for the case studies."""

from dataclasses import dataclass
from itertools import combinations
from math import isfinite
from numbers import Integral, Real as RealNumber

from z3 import (
    And,
    ForAll,
    Function,
    If,
    Implies,
    Int,
    IntSort,
    Not,
    Or,
    RealSort,
    RealVal,
    Solver,
    Sum,
    sat,
    unknown,
    unsat,
)


# The RV case-study encoding used two padded reference ticks per source tick.
# Keep that finite abstraction here so this remains the same full-flow
# baseline. Its completeness is checked against the archived continuous oracle.
_PAD = 2
_SOLVER_TIMEOUT_MILLISECONDS = None


@dataclass(frozen=True)
class _FullFlowEncoding:
    solver: Solver
    signal_functions: tuple[tuple[object, ...], ...]
    clocks: tuple[object, ...]
    variations: tuple[tuple[int, ...], ...]
    duration: int
    lower_tick: int
    upper_tick: int
    pad: int


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


def _validate_threshold(threshold, name, nonnegative=False):
    if (
        isinstance(threshold, bool)
        or not isinstance(threshold, RealNumber)
        or not isfinite(float(threshold))
        or (nonnegative and threshold < 0)
    ):
        qualifier = "finite nonnegative" if nonnegative else "finite"
        raise ValueError(f"{name} threshold must be a {qualifier} real number")


def _normalize_signals(signals, value_count):
    """Validate dense integer-grid PWC traces and discard terminal values."""
    if not 2 <= len(signals) <= 4:
        raise ValueError("case studies require two, three, or four agents")

    normalized = []
    duration = None
    for agent, data in enumerate(signals):
        if len(data) < 2:
            raise ValueError(
                f"data_{agent} must contain a sample and terminal marker"
            )

        rows = []
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
                or timestamp != index
            ):
                raise ValueError(
                    "case-study timestamps must be the dense integer grid "
                    "0,...,d"
                )

            values = []
            for value in row[1:value_count + 1]:
                if (
                    isinstance(value, bool)
                    or not isinstance(value, RealNumber)
                    or not isfinite(float(value))
                ):
                    raise ValueError(
                        "case-study values must be finite real numbers"
                    )
                values.append(value)
            rows.append(tuple(values))

        signal_duration = len(rows) - 1
        if duration is None:
            duration = signal_duration
        elif duration != signal_duration:
            raise ValueError(
                "all case-study traces must have the same duration"
            )
        normalized.append(tuple(rows[:-1]))

    return tuple(normalized), duration


def _variation_points(data):
    points = [0]
    for local_time in range(1, len(data)):
        if data[local_time] != data[local_time - 1]:
            points.append(local_time)
    points.append(len(data))
    return tuple(points)


def _allowed_local_times(tick, duration, padded_epsilon, pad):
    """Project the open reference-clock skew window to PWC indices."""
    return tuple(
        sorted(
            {
                min(
                    duration - 1,
                    max(
                        0,
                        (tick - padded_epsilon + offset) // pad,
                    ),
                )
                for offset in range(1, 2 * padded_epsilon)
            }
        )
    )


def _add_signal_functions(solver, rows, value_count):
    functions = []
    for agent, data in enumerate(rows):
        agent_functions = []
        for column in range(value_count):
            signal = Function(
                f"case_signal_{agent}_{column}", IntSort(), RealSort()
            )
            solver.add(
                *(
                    signal(local_time)
                    == RealVal(str(data[local_time][column]))
                    for local_time in range(len(data))
                )
            )
            agent_functions.append(signal)
        functions.append(tuple(agent_functions))
    return tuple(functions)


def _add_clock(
    solver,
    name,
    duration,
    lower_tick,
    upper_tick,
    epsilon,
    pad,
    variations,
):
    """Add one RV-style padded forward clock over the complete trace."""
    clock = Function(name, IntSort(), IntSort())
    padded_epsilon = pad * epsilon
    allowed = {
        tick: _allowed_local_times(
            tick, duration, padded_epsilon, pad
        )
        for tick in range(lower_tick, upper_tick + 1)
    }

    # This is the missing reference-clock constraint in the RV case scripts.
    # The window is open at both ends, reflecting strict continuous skew.
    solver.add(
        *(
            Or(*(clock(tick) == value for value in allowed[tick]))
            for tick in range(lower_tick, upper_tick + 1)
        )
    )
    solver.add(clock(lower_tick) == 0, clock(upper_tick) == duration - 1)

    # The PWC projection must visit every source interval exactly in order.
    for local_time in range(duration):
        candidate_ticks = (
            tick
            for tick in range(lower_tick, upper_tick + 1)
            if local_time in allowed[tick]
        )
        solver.add(
            Or(*(clock(tick) == local_time for tick in candidate_ticks))
        )

    # Retain the RV full-flow ordering workload, expressed without tautological
    # implications for pairs whose order is already known in Python.
    solver.add(
        *(
            clock(left) <= clock(right)
            for left in range(lower_tick, upper_tick + 1)
            for right in range(left + 1, upper_tick + 1)
        )
    )

    # Preserve each maximal constant piece, as in the RV adaptation.
    for index in range(len(variations) - 1):
        start = variations[index]
        end = variations[index + 1]
        solver.add(
            Or(
                *(
                    And(start <= clock(tick), clock(tick) < end)
                    for tick in range(lower_tick, upper_tick + 1)
                )
            )
        )
    return clock


def _edge_occurrence(clock, local_time, lower_tick, upper_tick):
    """Return the padded rank at which a PWC edge has occurred."""
    return Sum(
        *(
            If(clock(tick) < local_time, 1, 0)
            for tick in range(lower_tick, upper_tick + 1)
        )
    )


def _add_happened_before_constraints(encoding, epsilon):
    occurrences = []
    for clock, variations in zip(encoding.clocks, encoding.variations):
        occurrences.append(
            {
                local_time: _edge_occurrence(
                    clock,
                    local_time,
                    encoding.lower_tick,
                    encoding.upper_tick,
                )
                for local_time in variations[1:-1]
            }
        )

    for left, right in combinations(range(len(encoding.clocks)), 2):
        for left_time, left_occurrence in occurrences[left].items():
            for right_time, right_occurrence in occurrences[right].items():
                if left_time + epsilon <= right_time:
                    encoding.solver.add(left_occurrence < right_occurrence)
                elif right_time + epsilon <= left_time:
                    encoding.solver.add(right_occurrence < left_occurrence)


def _build_full_flow(epsilon, signals, value_count, pad=_PAD):
    _validate_epsilon(epsilon)
    rows, duration = _normalize_signals(signals, value_count)
    solver = _new_solver()
    lower_tick = 0
    upper_tick = duration * pad - 1
    variations = tuple(_variation_points(data) for data in rows)
    signal_functions = _add_signal_functions(solver, rows, value_count)
    clocks = tuple(
        _add_clock(
            solver,
            f"case_clock_{agent}",
            duration,
            lower_tick,
            upper_tick,
            epsilon,
            pad,
            variations[agent],
        )
        for agent in range(len(rows))
    )

    # Integer PWC indices can differ by exactly eps while the underlying real
    # readings remain strictly closer than eps. Equality is therefore allowed
    # here; strictness at the boundary is recovered by the open reference
    # domains and happened-before occurrence ordering.
    for left, right in combinations(clocks, 2):
        solver.add(
            *(
                And(
                    left(tick) - right(tick) <= epsilon,
                    right(tick) - left(tick) <= epsilon,
                )
                for tick in range(lower_tick, upper_tick + 1)
            )
        )

    encoding = _FullFlowEncoding(
        solver,
        signal_functions,
        clocks,
        variations,
        duration,
        lower_tick,
        upper_tick,
        pad,
    )
    _add_happened_before_constraints(encoding, epsilon)
    return encoding


def _values_at(encoding, tick):
    return tuple(
        tuple(signal(clock(tick)) for signal in agent_signals)
        for agent_signals, clock in zip(
            encoding.signal_functions, encoding.clocks
        )
    )


def _water_flow(encoding):
    """Materialize the RV aggregate flow over every padded tick."""
    flow = Function("case_water_flow", IntSort(), RealSort())
    encoding.solver.add(
        *(
            flow(tick)
            == Sum(
                *(agent[0] for agent in _values_at(encoding, tick))
            )
            for tick in range(
                encoding.lower_tick, encoding.upper_tick + 1
            )
        )
    )
    return flow


def _minimum(expressions):
    result = expressions[0]
    for expression in expressions[1:]:
        result = If(result <= expression, result, expression)
    return result


def _separation_flow(encoding):
    """Materialize the minimum pairwise squared distance at each tick."""
    flow = Function("case_separation_flow", IntSort(), RealSort())
    for tick in range(encoding.lower_tick, encoding.upper_tick + 1):
        values = _values_at(encoding, tick)
        squared_distances = []
        for left, right in combinations(range(len(values)), 2):
            squared_distances.append(
                Sum(
                    *(
                        (values[left][axis] - values[right][axis])
                        * (values[left][axis] - values[right][axis])
                        for axis in range(3)
                    )
                )
            )
        encoding.solver.add(
            flow(tick) == _minimum(squared_distances)
        )
    return flow


def _check_unsat(solver):
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


def _prove_positive(encoding, property_at):
    witness = Int("case_bad_tick")
    encoding.solver.add(
        witness >= encoding.lower_tick,
        witness <= encoding.upper_tick,
        Not(property_at(witness)),
    )
    return _check_unsat(encoding.solver)


def _prove_negative(encoding, property_at):
    tick = Int("case_safe_tick")
    encoding.solver.add(
        ForAll(
            tick,
            Implies(
                And(
                    tick >= encoding.lower_tick,
                    tick <= encoding.upper_tick,
                ),
                property_at(tick),
            ),
        )
    )
    return _check_unsat(encoding.solver)


def prove_water_tank(epsilon, signals, threshold=10):
    """Prove G(sum(signals) > threshold) for every admissible flow."""
    _validate_threshold(threshold, "water-tank")
    encoding = _build_full_flow(epsilon, signals, 1)
    flow = _water_flow(encoding)
    threshold_value = RealVal(str(threshold))
    return _prove_positive(
        encoding, lambda tick: flow(tick) > threshold_value
    )


def prove_water_tank_negative(epsilon, signals, threshold=10):
    """Prove every admissible flow violates G(sum > threshold)."""
    _validate_threshold(threshold, "water-tank")
    encoding = _build_full_flow(epsilon, signals, 1)
    flow = _water_flow(encoding)
    threshold_value = RealVal(str(threshold))
    return _prove_negative(
        encoding, lambda tick: flow(tick) > threshold_value
    )


def prove_mutual_separation(epsilon, signals, threshold=0):
    """Prove every agent pair remains farther apart than threshold."""
    _validate_threshold(threshold, "separation", nonnegative=True)
    encoding = _build_full_flow(epsilon, signals, 3)
    flow = _separation_flow(encoding)
    threshold_value = RealVal(str(threshold))
    squared_threshold = threshold_value * threshold_value
    return _prove_positive(
        encoding, lambda tick: flow(tick) > squared_threshold
    )


def prove_mutual_separation_negative(epsilon, signals, threshold=0):
    """Prove every admissible flow contains a separation violation."""
    _validate_threshold(threshold, "separation", nonnegative=True)
    encoding = _build_full_flow(epsilon, signals, 3)
    flow = _separation_flow(encoding)
    threshold_value = RealVal(str(threshold))
    squared_threshold = threshold_value * threshold_value
    return _prove_negative(
        encoding, lambda tick: flow(tick) > squared_threshold
    )
