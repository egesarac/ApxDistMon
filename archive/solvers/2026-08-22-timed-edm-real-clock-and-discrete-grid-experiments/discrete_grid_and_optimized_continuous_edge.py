"""Archived discrete-grid and optimized continuous-edge timed encodings.

The *_timed entry points solve a padded discrete-grid clock model. That is a
different semantics from the paper's continuous-time retimings and must not
be used as the paper baseline. The *_timed_continuous entry points are the
experimental optimized continuous edge-based encoding.

This snapshot also contains the untimed helpers that shared the module when
the alternatives were maintained. No entry point here is used by the current
benchmark runner.
"""

from dataclasses import dataclass
from numbers import Integral, Real as RealNumber

from z3 import (
    And,
    BoolVal,
    ForAll,
    Function,
    Implies,
    If,
    Int,
    IntSort,
    Not,
    Or,
    Real,
    RealVal,
    Solver,
    Sum,
    sat,
    set_option,
    unknown,
    unsat,
)


set_option(rational_to_decimal=True)
set_option(
    max_args=100000000,
    max_lines=10000000000,
    max_depth=100000000,
    max_visited=10000000,
)


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


@dataclass(frozen=True)
class _DiscreteEncoding:
    solver: Solver
    sig0: object
    sig1: object
    c0: object
    c1: object
    lower_tick: int
    upper_tick: int


@dataclass(frozen=True)
class _TimedDiscreteEncoding:
    solver: Solver
    clock0: tuple[object, ...]
    clock1: tuple[object, ...]
    scale: int
    horizon: int


@dataclass(frozen=True)
class _Edge:
    agent: int
    index: int
    local_time: int
    value: bool


def _validate_epsilon(eps):
    if isinstance(eps, bool) or not isinstance(eps, Integral) or eps <= 0:
        raise ValueError("eps must be a positive integer")


def _check_solver(solver):
    result = solver.check()
    if result == unknown:
        reason = solver.reason_unknown()
        if "timeout" in reason.lower():
            raise TimeoutError(f"Z3 timed out: {reason}")
        raise RuntimeError(f"Z3 returned unknown: {reason}")
    if result not in (sat, unsat):
        raise RuntimeError(f"unexpected Z3 result: {result}")
    return result


def _variation_points(data):
    points = {0}
    for index in range(1, len(data)):
        if data[index][1] != data[index - 1][1]:
            points.add(int(data[index][0]))
    points.add(int(data[-1][0]))
    return points


def _add_discrete_signal(
    solver,
    name,
    data,
    variation_points,
    segment_lower_bound,
    segment_upper_bound,
    segment_start,
    signal_duration,
    pad,
    entry_found,
):
    signal = Function(name, IntSort(), IntSort())
    timestamps = []
    segment_variation = []

    for timestamp in range(segment_lower_bound, segment_upper_bound + 1):
        if 0 <= timestamp < len(data) - 1:
            solver.add(signal(timestamp) == data[timestamp][1])
            for offset in range(pad):
                timestamps.append(timestamp * pad + offset)

            if timestamp == segment_lower_bound or (
                data[timestamp][0] in variation_points
                and data[timestamp][0] > segment_lower_bound
            ):
                segment_variation.append(int(data[timestamp][0]))

            if timestamp > segment_start:
                entry_found = True
            if timestamp == signal_duration:
                entry_found = False

    if (
        segment_variation
        and segment_variation[-1] < segment_upper_bound
    ):
        segment_variation.append(segment_upper_bound)

    return signal, timestamps, segment_variation, entry_found


def _add_discrete_clock(
    solver,
    name,
    data,
    timestamps,
    segment_lower_bound,
    segment_upper_bound,
    pad,
    padded_epsilon,
    enumerate_domain,
    inverse_offsets,
):
    clock = Function(name, IntSort(), IntSort())
    first_tick = timestamps[0]
    last_tick = timestamps[-1]

    if enumerate_domain:
        solver.add(
            And(
                [
                    Or(
                        [
                            clock(tick)
                            == min(
                                last_tick // pad,
                                max(
                                    0,
                                    (tick - padded_epsilon + offset) // pad,
                                ),
                            )
                            for offset in range(1, 2 * int(padded_epsilon))
                        ]
                    )
                    for tick in range(first_tick, last_tick + 1)
                ]
            )
        )
    else:
        solver.add(
            And(
                [
                    And(clock(tick) >= 0, clock(tick) < data[-1][0])
                    for tick in range(first_tick, last_tick + 1)
                ]
            )
        )

    solver.add(
        And(
            [
                Or(
                    [
                        clock(
                            min(
                                last_tick,
                                max(
                                    0,
                                    local_time * pad
                                    - padded_epsilon
                                    + offset,
                                ),
                            )
                        )
                        == local_time
                        for offset in inverse_offsets
                    ]
                )
                for local_time in range(
                    segment_lower_bound, segment_upper_bound
                )
            ]
        )
    )

    if enumerate_domain:
        solver.add(
            And(
                clock(first_tick) >= segment_lower_bound,
                clock(last_tick) < segment_upper_bound,
            )
        )

    return clock


def _add_variation_constraints(
    solver, clock, segment_variation, first_tick, last_tick
):
    solver.add(
        And(
            [
                Or(
                    [
                        And(
                            segment_variation[index] <= clock(tick),
                            clock(tick) < segment_variation[index + 1],
                        )
                        for tick in range(first_tick, last_tick + 1)
                    ]
                )
                for index in range(len(segment_variation) - 1)
            ]
        )
    )


def _build_discrete_encoding(
    eps,
    data_0,
    data_1,
    variation_0,
    variation_1,
    signal_duration,
    segment_duration,
    segment_index,
    enumerate_domain,
):
    pad = 2
    padded_epsilon = pad * eps
    segment_lower_bound = max(
        0, int(segment_index * segment_duration - eps)
    )
    segment_upper_bound = min(
        signal_duration,
        int((segment_index + 1) * segment_duration + eps),
    )

    solver = _new_solver()
    entry_found = False
    sig0, timestamps_0, segment_variation_0, entry_found = (
        _add_discrete_signal(
            solver,
            "sig0",
            data_0,
            variation_0,
            segment_lower_bound,
            segment_upper_bound,
            segment_index * segment_duration,
            signal_duration,
            pad,
            entry_found,
        )
    )
    sig1, timestamps_1, segment_variation_1, entry_found = (
        _add_discrete_signal(
            solver,
            "sig1",
            data_1,
            variation_1,
            segment_lower_bound,
            segment_upper_bound,
            segment_index * segment_duration,
            signal_duration,
            pad,
            entry_found,
        )
    )

    c0 = _add_discrete_clock(
        solver,
        "c0",
        data_0,
        timestamps_0,
        segment_lower_bound,
        segment_upper_bound,
        pad,
        padded_epsilon,
        enumerate_domain,
        range(1, 2 * int(padded_epsilon) + pad),
    )
    c1 = _add_discrete_clock(
        solver,
        "c1",
        data_1,
        timestamps_1,
        segment_lower_bound,
        segment_upper_bound,
        pad,
        padded_epsilon,
        enumerate_domain,
        range(2 * int(padded_epsilon) + pad),
    )

    first_tick = timestamps_0[0]
    last_tick = timestamps_0[-1]
    _add_variation_constraints(
        solver, c0, segment_variation_0, first_tick, last_tick
    )
    _add_variation_constraints(
        solver, c1, segment_variation_1, first_tick, last_tick
    )

    solver.add(
        And(
            [
                And(c0(tick) - c1(tick) <= eps, c0(tick) - c1(tick) >= -eps)
                for tick in range(first_tick, last_tick + 1)
            ]
        )
    )
    solver.add(
        And(
            [
                Implies(
                    And(
                        Or(
                            [
                                And(
                                    c0(tick) == segment_variation_0[index],
                                    c0(tick) != c0(tick - 1),
                                )
                                for index in range(
                                    1, len(segment_variation_0) - 1
                                )
                            ]
                        ),
                        Or(
                            [
                                And(
                                    c1(tick) == segment_variation_1[index],
                                    c1(tick) != c1(tick - 1),
                                )
                                for index in range(
                                    1, len(segment_variation_1) - 1
                                )
                            ]
                        ),
                    ),
                    And(
                        c0(tick) - c1(tick) < eps,
                        c0(tick) - c1(tick) > -eps,
                    ),
                )
                for tick in range(timestamps_1[0], timestamps_1[-1] + 1)
            ]
        )
    )
    solver.add(
        And(
            [
                And(
                    [
                        Implies(
                            left <= right,
                            And(c0(left) <= c0(right), c1(left) <= c1(right)),
                        )
                        for right in range(first_tick, last_tick + 1)
                    ]
                )
                for left in range(first_tick, last_tick + 1)
            ]
        )
    )

    return (
        _DiscreteEncoding(
            solver, sig0, sig1, c0, c1, first_tick, last_tick
        ),
        entry_found,
    )


def _consistent_cut_flow(encoding):
    flow = Function("c_flow", IntSort(), IntSort())
    encoding.solver.add(
        And(
            [
                flow(tick)
                == encoding.sig0(encoding.c0(tick))
                + encoding.sig1(encoding.c1(tick))
                for tick in range(
                    encoding.lower_tick, encoding.upper_tick + 1
                )
            ]
        )
    )
    return flow


def _add_untimed_query(encoding, query):
    solver = encoding.solver
    lower_tick = encoding.lower_tick
    upper_tick = encoding.upper_tick

    if query == "not_response":
        # Kept for parity with the original encoding, where this definition
        # was present even though the temporal query does not reference it.
        _consistent_cut_flow(encoding)
        u = Int("u")
        v = Int("v")
        solver.add(And(v >= u, v <= upper_tick))
        solver.add(And(u >= lower_tick, u <= upper_tick))
        solver.add(
            ForAll(
                u,
                Implies(
                    And(
                        u >= lower_tick,
                        u <= upper_tick,
                        encoding.sig0(encoding.c0(u)) == 1,
                    ),
                    And(
                        v >= u,
                        v <= upper_tick,
                        encoding.sig1(encoding.c1(v)) == 1,
                    ),
                ),
            )
        )
        return

    if query == "response":
        u = Int("u")
        v = Int("v")
        solver.add(And(v >= u, v <= upper_tick))
        solver.add(And(u >= lower_tick, u <= upper_tick))
        solver.add(
            And(
                u >= lower_tick,
                u <= upper_tick,
                encoding.sig0(encoding.c0(u)) == 1,
                ForAll(
                    v,
                    Implies(
                        And(v >= u, v <= upper_tick),
                        encoding.sig1(encoding.c1(v)) == 0,
                    ),
                ),
            )
        )
        return

    if query in {
        "always_conjunction",
        "always_disjunction",
        "eventually_conjunction",
        "eventually_disjunction",
    }:
        flow = _consistent_cut_flow(encoding)
        v = Int("v")
        bounds = And(v >= lower_tick, v <= upper_tick)
        solver.add(bounds)
        if query.endswith("conjunction"):
            state_formula = flow(v) < 2
        else:
            state_formula = flow(v) == 0
        if query.startswith("eventually"):
            solver.add(ForAll(v, Implies(bounds, state_formula)))
        else:
            solver.add(Implies(bounds, state_formula))
        return

    if query == "until":
        _consistent_cut_flow(encoding)
        v = Int("v")
        solver.add(And(v >= lower_tick, v <= upper_tick))
        u = Int("u")
        solver.add(And(u >= lower_tick, u <= v))
        solver.add(
            ForAll(
                v,
                Implies(
                    And(v >= lower_tick, v <= upper_tick),
                    Or(
                        encoding.sig1(encoding.c1(v)) == 0,
                        And(
                            u >= lower_tick,
                            u < v,
                            encoding.sig0(encoding.c0(u)) == 0,
                        ),
                    ),
                ),
            )
        )
        return

    if query == "not_until":
        _consistent_cut_flow(encoding)
        v = Int("v")
        solver.add(And(v >= lower_tick, v <= upper_tick))
        u = Int("u")
        solver.add(And(u >= lower_tick, u <= v))
        solver.add(
            And(
                v >= lower_tick,
                v <= upper_tick,
                encoding.sig1(encoding.c1(v)) == 1,
                ForAll(
                    u,
                    Implies(
                        And(u >= lower_tick, u < v),
                        encoding.sig0(encoding.c0(u)) == 1,
                    ),
                ),
            )
        )
        return

    raise ValueError(f"unsupported untimed query: {query}")


def _prove_untimed(eps, segCount, data_0, data_1, query):
    _validate_epsilon(eps)

    signal_duration = len(data_0) - 1
    segment_duration = signal_duration / segCount
    t0 = data_0[0][0]
    t1 = data_0[1][0]
    if signal_duration / segCount < t1:
        segCount = signal_duration / t1
    if t0 != 0:
        return None

    variation_0 = _variation_points(data_0)
    variation_1 = _variation_points(data_1)
    enumerate_domain = query in {
        "not_response",
        "response",
        "until",
        "not_until",
    }
    stop_on_sat = query in {"always_conjunction", "always_disjunction"}
    stop_on_unsat = query in {
        "eventually_conjunction",
        "eventually_disjunction",
    }

    segment_index = 0
    entry_found = True
    flag = True
    while entry_found:
        if segment_index == segCount:
            return flag

        encoding, entry_found = _build_discrete_encoding(
            eps,
            data_0,
            data_1,
            variation_0,
            variation_1,
            signal_duration,
            segment_duration,
            segment_index,
            enumerate_domain,
        )
        segment_index += 1
        _add_untimed_query(encoding, query)

        result = _check_solver(encoding.solver)
        if result == sat:
            flag = False
            if stop_on_sat:
                return flag
        elif segment_index <= segCount:
            flag = True
            if stop_on_unsat:
                return flag

    return flag


def prog_not_always_implies_eventually(eps, segCount, data_0, data_1):
    return _prove_untimed(
        eps, segCount, data_0, data_1, "not_response"
    )


def prog_always_implies_eventually(eps, segCount, data_0, data_1):
    return _prove_untimed(eps, segCount, data_0, data_1, "response")


def prog_always_conjunction(eps, segCount, data_0, data_1):
    return _prove_untimed(
        eps, segCount, data_0, data_1, "always_conjunction"
    )


def prog_always_disjunction(eps, segCount, data_0, data_1):
    return _prove_untimed(
        eps, segCount, data_0, data_1, "always_disjunction"
    )


def prog_eventually_conjunction(eps, segCount, data_0, data_1):
    return _prove_untimed(
        eps, segCount, data_0, data_1, "eventually_conjunction"
    )


def prog_eventually_disjunction(eps, segCount, data_0, data_1):
    return _prove_untimed(
        eps, segCount, data_0, data_1, "eventually_disjunction"
    )


def prog_until(eps, segCount, data_0, data_1):
    return _prove_untimed(eps, segCount, data_0, data_1, "until")


def prog_not_until(eps, segCount, data_0, data_1):
    return _prove_untimed(eps, segCount, data_0, data_1, "not_until")


def _validate_bounds(a, b):
    if (
        isinstance(a, bool)
        or not isinstance(a, Integral)
        or isinstance(b, bool)
        or not isinstance(b, Integral)
    ):
        raise ValueError("response bounds must be integers")
    if a < 0 or a >= b:
        raise ValueError("response bounds must satisfy 0 <= a < b")


def _validate_signal(data, name):
    if len(data) < 2:
        raise ValueError(f"{name} must contain a sample and a terminal marker")

    for index, row in enumerate(data):
        if len(row) < 2:
            raise ValueError(f"{name}[{index}] must contain time and value")
        timestamp, value = row[:2]
        if (
            isinstance(timestamp, bool)
            or not isinstance(timestamp, RealNumber)
            or timestamp != index
        ):
            raise ValueError(
                f"{name} timestamps must be the integer grid 0,...,d"
            )
        if value not in (0, 1):
            raise ValueError(f"{name}[{index}] must have Boolean value 0 or 1")

    if data[-1][1] != data[-2][1]:
        raise ValueError(
            f"{name} must repeat its final value at the terminal marker"
        )

    return len(data) - 1


def _validate_timed_inputs(eps, data_0, data_1, a, b):
    _validate_epsilon(eps)
    _validate_bounds(a, b)
    duration_0 = _validate_signal(data_0, "data_0")
    duration_1 = _validate_signal(data_1, "data_1")
    if duration_0 != duration_1:
        raise ValueError("data_0 and data_1 must have the same duration")
    return duration_0


def _edges(data, agent):
    return [
        _Edge(agent, index, index, data[index][1] == 1)
        for index in range(1, len(data) - 1)
        if data[index][1] != data[index - 1][1]
    ]


def _or(expressions):
    return Or(*expressions) if expressions else BoolVal(False)


def _and(expressions):
    return And(*expressions) if expressions else BoolVal(True)


def _timed_grid_scale(eps):
    """Return the number of global ticks per local-time unit."""
    # This is part of the discrete semantics and must not be capped by the
    # signal duration, as it was in the archived implementation.
    return 2 * eps


def _padded_cell_bounds(tick, scale, eps, horizon):
    """Return local cells intersecting the strict epsilon neighbourhood."""
    padded_epsilon = scale * eps
    fine_lower = max(0, tick - padded_epsilon + 1)
    fine_upper = min(horizon - 1, tick + padded_epsilon - 1)
    return fine_lower // scale, fine_upper // scale


def _add_timed_discrete_clock(
    solver, name, duration, eps, scale, horizon
):
    clock = tuple(Int(f"{name}_{tick}") for tick in range(horizon))

    # Starting at zero, ending in the final cell, and advancing by at most
    # one cell makes the clock monotone and surjective without the original
    # quadratic order constraints or large inverse-domain disjunctions.
    solver.add(clock[0] == 0, clock[-1] == duration - 1)
    for tick, reading in enumerate(clock):
        cell_lower, cell_upper = _padded_cell_bounds(
            tick, scale, eps, horizon
        )
        solver.add(cell_lower <= reading, reading <= cell_upper)

    for left, right in zip(clock, clock[1:]):
        solver.add(left <= right, right <= left + 1)

    return clock


def _timed_variation_indices(data):
    return tuple(
        index
        for index in range(1, len(data) - 1)
        if data[index][1] != data[index - 1][1]
    )


def _add_timed_edge_occurrences(solver, name, clock, variations):
    occurrences = {}
    for local_time in variations:
        occurrence = Int(f"{name}_{local_time}")
        # A monotone, surjective unit-step clock first reaches local_time
        # after exactly this many readings below local_time.
        solver.add(
            occurrence
            == Sum(
                *[
                    If(reading < local_time, 1, 0)
                    for reading in clock
                ]
            )
        )
        occurrences[local_time] = occurrence
    return occurrences


def _add_timed_happened_before(
    solver, eps, occurrences_0, occurrences_1
):
    # Events whose local timestamps differ by at least epsilon cannot be
    # concurrent. Events closer than epsilon retain both possible orders.
    for local_time_0, occurrence_0 in occurrences_0.items():
        for local_time_1, occurrence_1 in occurrences_1.items():
            if local_time_0 + eps <= local_time_1:
                solver.add(occurrence_0 < occurrence_1)
            elif local_time_1 + eps <= local_time_0:
                solver.add(occurrence_1 < occurrence_0)


def _build_timed_discrete_encoding(eps, data_0, data_1, duration):
    scale = _timed_grid_scale(eps)
    horizon = scale * duration
    solver = _new_solver()
    clock0 = _add_timed_discrete_clock(
        solver, "timed_c0", duration, eps, scale, horizon
    )
    clock1 = _add_timed_discrete_clock(
        solver, "timed_c1", duration, eps, scale, horizon
    )

    # Cell indices use inclusive epsilon skew. The fine-tick neighbourhood
    # above remains strict; exact-boundary event order is handled by HB below.
    for reading_0, reading_1 in zip(clock0, clock1):
        solver.add(
            -eps <= reading_0 - reading_1,
            reading_0 - reading_1 <= eps,
        )

    occurrences_0 = _add_timed_edge_occurrences(
        solver,
        "timed_occurrence_0",
        clock0,
        _timed_variation_indices(data_0),
    )
    occurrences_1 = _add_timed_edge_occurrences(
        solver,
        "timed_occurrence_1",
        clock1,
        _timed_variation_indices(data_1),
    )
    _add_timed_happened_before(
        solver, eps, occurrences_0, occurrences_1
    )

    return _TimedDiscreteEncoding(
        solver, clock0, clock1, scale, horizon
    )


def _true_runs(data):
    runs = []
    start = None
    for index, (_, value) in enumerate(data[:-1]):
        if value == 1 and start is None:
            start = index
        elif value == 0 and start is not None:
            runs.append((start, index))
            start = None

    if start is not None:
        runs.append((start, len(data) - 1))
    return tuple(runs)


def _timed_signal_value(clock_reading, true_runs):
    return _or(
        [
            And(start <= clock_reading, clock_reading < end)
            for start, end in true_runs
        ]
    )


def _timed_discrete_response_formula(encoding, data_0, data_1, a, b):
    p_runs = _true_runs(data_0)
    q_runs = _true_runs(data_1)
    p_values = tuple(
        _timed_signal_value(reading, p_runs)
        for reading in encoding.clock0
    )
    q_values = tuple(
        _timed_signal_value(reading, q_runs)
        for reading in encoding.clock1
    )

    # Prefix counts share each q-value across all response windows, avoiding
    # the original quantified witness and an explicit quadratic expansion.
    q_prefix = tuple(
        Int(f"timed_q_prefix_{tick}")
        for tick in range(encoding.horizon + 1)
    )
    encoding.solver.add(q_prefix[0] == 0)
    for tick, q_value in enumerate(q_values):
        encoding.solver.add(
            q_prefix[tick + 1]
            == q_prefix[tick] + If(q_value, 1, 0)
        )

    lower_offset = a * encoding.scale
    upper_offset = b * encoding.scale
    obligations = []
    for tick, p_value in enumerate(p_values):
        window_lower = min(encoding.horizon, tick + lower_offset)
        window_upper = min(encoding.horizon, tick + upper_offset)
        has_witness = (
            q_prefix[window_upper] > q_prefix[window_lower]
            if window_lower < window_upper
            else BoolVal(False)
        )
        obligations.append(Implies(p_value, has_witness))

    return _and(obligations)


def _prove_bounded_response_discrete(
    eps, data_0, data_1, a, b, prove_negation=False
):
    duration = _validate_timed_inputs(eps, data_0, data_1, a, b)
    encoding = _build_timed_discrete_encoding(
        eps, data_0, data_1, duration
    )
    response = _timed_discrete_response_formula(
        encoding, data_0, data_1, a, b
    )

    encoding.solver.add(response if prove_negation else Not(response))
    return _check_solver(encoding.solver) == unsat


def _continuous_clock_model(eps, data_0, data_1, duration):
    edges_by_agent = (_edges(data_0, 0), _edges(data_1, 1))
    edges = edges_by_agent[0] + edges_by_agent[1]
    solver = _new_solver()
    duration_value = RealVal(duration)
    epsilon_value = RealVal(eps)

    occurrence = {
        edge: Real(f"occurrence_{edge.agent}_{edge.index}") for edge in edges
    }
    reading = {
        (agent, edge): Real(
            f"reading_{agent}_{edge.agent}_{edge.index}"
        )
        for agent in (0, 1)
        for edge in edges
    }

    for edge in edges:
        physical_time = occurrence[edge]
        solver.add(0 < physical_time, physical_time < duration_value)
        solver.add(reading[edge.agent, edge] == edge.local_time)

        for agent in (0, 1):
            local_time = reading[agent, edge]
            solver.add(0 < local_time, local_time < duration_value)
            solver.add(local_time - physical_time < epsilon_value)
            solver.add(physical_time - local_time < epsilon_value)

        solver.add(
            reading[0, edge] - reading[1, edge] < epsilon_value
        )
        solver.add(
            reading[1, edge] - reading[0, edge] < epsilon_value
        )

    for index, left in enumerate(edges):
        for right in edges[index + 1 :]:
            left_time = occurrence[left]
            right_time = occurrence[right]
            solver.add(
                Or(
                    And(
                        left_time < right_time,
                        reading[0, left] < reading[0, right],
                        reading[1, left] < reading[1, right],
                    ),
                    And(
                        left_time == right_time,
                        reading[0, left] == reading[0, right],
                        reading[1, left] == reading[1, right],
                    ),
                    And(
                        left_time > right_time,
                        reading[0, left] > reading[0, right],
                        reading[1, left] > reading[1, right],
                    ),
                )
            )

    return solver, edges_by_agent, occurrence


def _constant_intervals(data, agent_edges, occurrence, wanted, duration):
    bounds = (
        [RealVal(0)]
        + [occurrence[edge] for edge in agent_edges]
        + [RealVal(duration)]
    )
    values = [data[0][1] == 1] + [edge.value for edge in agent_edges]
    return [
        (bounds[index], bounds[index + 1])
        for index, value in enumerate(values)
        if value == wanted
    ]


def _response_formula(
    data_0, data_1, edges_by_agent, occurrence, duration, a, b
):
    p_intervals = _constant_intervals(
        data_0, edges_by_agent[0], occurrence, True, duration
    )
    q_intervals = _constant_intervals(
        data_1, edges_by_agent[1], occurrence, True, duration
    )
    q_windows = [(start - b, end - a) for start, end in q_intervals]

    p_coverage = []
    for p_start, p_end in p_intervals:
        # Whether a connected chain ending at the current q-window
        # contains the start of this p-interval.
        reaches_start = BoolVal(False)
        covers_interval = []
        previous_right = None

        for left, right in q_windows:
            starts_at_current = And(left < p_start, p_start < right)
            if previous_right is None:
                reaches_start = starts_at_current
            else:
                reaches_start = Or(
                    starts_at_current,
                    And(reaches_start, left < previous_right),
                )

            covers_interval.append(
                And(reaches_start, right >= p_end)
            )
            previous_right = right

        p_coverage.append(_or(covers_interval))

    return _and(p_coverage)


def _prove_bounded_response_continuous(
    eps, data_0, data_1, a, b, prove_negation=False
):
    duration = _validate_timed_inputs(eps, data_0, data_1, a, b)
    solver, edges_by_agent, occurrence = _continuous_clock_model(
        eps, data_0, data_1, duration
    )
    response = _response_formula(
        data_0, data_1, edges_by_agent, occurrence, duration, a, b
    )

    solver.add(response if prove_negation else Not(response))
    result = _check_solver(solver)
    return result == unsat


def prog_not_always_implies_eventually_timed_continuous(
    eps, segCount, data_0, data_1, a, b
):
    """Continuous-time reference: prove every trace violates response.

    segCount is retained for compatibility; timed formulas use the full trace.
    """
    return _prove_bounded_response_continuous(
        eps, data_0, data_1, a, b, prove_negation=True
    )


def prog_always_implies_eventually_timed_continuous(
    eps, segCount, data_0, data_1, a, b
):
    """Continuous-time reference: prove every trace satisfies response.

    segCount is retained for compatibility; timed formulas use the full trace.
    """
    return _prove_bounded_response_continuous(
        eps, data_0, data_1, a, b
    )


def prog_not_always_implies_eventually_timed(
    eps, segCount, data_0, data_1, a, b
):
    """Prove every full-horizon discrete retiming violates response.

    segCount is retained for compatibility; timed formulas use the full trace.
    """
    return _prove_bounded_response_discrete(
        eps, data_0, data_1, a, b, prove_negation=True
    )


def prog_always_implies_eventually_timed(
    eps, segCount, data_0, data_1, a, b
):
    """Prove every full-horizon discrete retiming satisfies response.

    segCount is retained for compatibility; timed formulas use the full trace.
    """
    return _prove_bounded_response_discrete(
        eps, data_0, data_1, a, b
    )


def preprocess(data, d):
    out = []
    for index in range(d):
        if data[index][1] > 0:
            out.append([data[index][0], 1.0])
        else:
            out.append([data[index][0], 0.0])
    out.append([float(d), out[d - 1][1]])
    return out


def negate(data):
    return [[row[0], 1.0 - row[1]] for row in data]
