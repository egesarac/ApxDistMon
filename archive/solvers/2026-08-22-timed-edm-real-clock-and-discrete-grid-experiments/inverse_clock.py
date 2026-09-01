"""Archived inverse-clock experiment for exact timed monitoring.

This exploratory encoding represents continuous clocks through real-valued
inverse-clock knots on the integer local-time grid. It is retained only for
historical comparison and is not used by the maintained benchmark runner.
"""

from dataclasses import dataclass
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
    Real,
    RealSort,
    RealVal,
    Solver,
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
class _TimedGridEncoding:
    solver: Solver
    sig0: object
    sig1: object
    inverse0: object
    inverse1: object
    duration: int


def _validate_epsilon(eps):
    if isinstance(eps, bool) or not isinstance(eps, Integral) or eps <= 0:
        raise ValueError("eps must be a positive integer")


def _check_solver(solver):
    result = solver.check()
    if result == unknown:
        raise RuntimeError(f"Z3 returned unknown: {solver.reason_unknown()}")
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

    solver = Solver()
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


def _validate_timed_inputs(eps, seg_count, data_0, data_1, a, b):
    _validate_epsilon(eps)
    _validate_bounds(a, b)
    if seg_count != 1:
        raise ValueError(
            "the offline exact timed encoding requires segCount == 1"
        )

    duration_0 = _validate_signal(data_0, "data_0")
    duration_1 = _validate_signal(data_1, "data_1")
    if duration_0 != duration_1:
        raise ValueError("data_0 and data_1 must have the same duration")
    return duration_0


def _add_timed_signal(solver, name, data, duration):
    signal = Function(name, IntSort(), IntSort())
    solver.add(
        And(
            [
                signal(local_time) == int(data[local_time][1])
                for local_time in range(duration + 1)
            ]
        )
    )
    return signal


def _add_inverse_clock(solver, name, duration):
    """Add real inverse-clock values on the integer local-time grid."""
    inverse = Function(name, IntSort(), RealSort())
    solver.add(inverse(0) == 0, inverse(duration) == duration)
    return inverse


def _build_timed_grid_encoding(eps, data_0, data_1, duration):
    """Encode an admissible clock family by inverse-clock grid knots."""
    solver = Solver()
    sig0 = _add_timed_signal(
        solver, "timed_sig0", data_0, duration
    )
    sig1 = _add_timed_signal(
        solver, "timed_sig1", data_1, duration
    )
    inverse0 = _add_inverse_clock(solver, "timed_inverse0", duration)
    inverse1 = _add_inverse_clock(solver, "timed_inverse1", duration)

    # Inverse-clock mappings are strictly ordered. This deliberately follows
    # the original baseline's all-pairs grid encoding rather than replacing it
    # by the equivalent adjacent constraints.
    solver.add(
        And(
            [
                And(
                    [
                        Implies(
                            left < right,
                            And(
                                inverse0(left) < inverse0(right),
                                inverse1(left) < inverse1(right),
                            ),
                        )
                        for right in range(duration + 1)
                    ]
                )
                for left in range(duration + 1)
            ]
        )
    )

    # The first clock is the inverse of the semantic reference clock.
    clocks = [
        [RealVal(local_time) for local_time in range(duration + 1)],
        [inverse0(local_time) for local_time in range(duration + 1)],
        [inverse1(local_time) for local_time in range(duration + 1)],
    ]

    # For integer eps, g_i(t) < g_j(t + eps) at every integer t is
    # equivalent to strict eps-skew for the piecewise-linear inverse clocks.
    if eps <= duration:
        for source in range(len(clocks)):
            for target in range(len(clocks)):
                if source == target:
                    continue
                solver.add(
                    And(
                        [
                            clocks[source][local_time]
                            < clocks[target][local_time + eps]
                            for local_time in range(duration - eps + 1)
                        ]
                    )
                )

    return _TimedGridEncoding(
        solver, sig0, sig1, inverse0, inverse1, duration
    )


def _z3_step_value(signal, inverse, physical_time, duration):
    """Evaluate a right-continuous Boolean signal at a real physical time."""
    value = signal(0)
    for local_time in range(1, duration):
        value = If(
            physical_time >= inverse(local_time),
            signal(local_time),
            value,
        )
    return value


def _add_timed_response_query(encoding, a, b, satisfying_flow):
    """Add either a satisfying-flow or a violating-flow response query."""
    duration = RealVal(encoding.duration)
    lower = RealVal(a)
    upper = RealVal(b)
    u = Real("timed_response_u")

    p_holds = (
        _z3_step_value(
            encoding.sig0,
            encoding.inverse0,
            u,
            encoding.duration,
        )
        == 1
    )
    trigger = And(0 <= u, u < duration, p_holds)

    if satisfying_flow:
        # The free function is the Skolem witness for the existential
        # time in F_[a,b) q; unlike one shared free variable, it may choose
        # a different witness for every trigger time.
        witness = Function(
            "timed_response_witness", RealSort(), RealSort()
        )
        witness_time = witness(u)
        q_holds = (
            _z3_step_value(
                encoding.sig1,
                encoding.inverse1,
                witness_time,
                encoding.duration,
            )
            == 1
        )
        encoding.solver.add(
            ForAll(
                [u],
                Implies(
                    trigger,
                    And(
                        0 <= witness_time,
                        witness_time < duration,
                        u + lower <= witness_time,
                        witness_time < u + upper,
                        q_holds,
                    ),
                ),
            )
        )
        return

    # A violating flow has one trigger time for which q is false at every
    # point of [u+a,u+b) that still lies in the finite trace.
    v = Real("timed_response_v")
    q_holds = (
        _z3_step_value(
            encoding.sig1,
            encoding.inverse1,
            v,
            encoding.duration,
        )
        == 1
    )
    encoding.solver.add(trigger)
    encoding.solver.add(
        ForAll(
            [v],
            Implies(
                And(
                    0 <= v,
                    v < duration,
                    u + lower <= v,
                    v < u + upper,
                ),
                Not(q_holds),
            ),
        )
    )

def _prove_bounded_response(
    eps, seg_count, data_0, data_1, a, b, prove_negation=False
):
    duration = _validate_timed_inputs(
        eps, seg_count, data_0, data_1, a, b
    )
    encoding = _build_timed_grid_encoding(
        int(eps), data_0, data_1, duration
    )

    # For universal satisfaction, search for a violating clock family.
    # For universal violation, search for a satisfying clock family.
    _add_timed_response_query(
        encoding, a, b, satisfying_flow=prove_negation
    )
    result = _check_solver(encoding.solver)
    return result == unsat


def prog_not_always_implies_eventually_timed(
    eps, segCount, data_0, data_1, a, b
):
    """Return whether every admissible trace violates bounded response."""
    return _prove_bounded_response(
        eps,
        segCount,
        data_0,
        data_1,
        a,
        b,
        prove_negation=True,
    )


def prog_always_implies_eventually_timed(
    eps, segCount, data_0, data_1, a, b
):
    """Return whether every admissible trace satisfies bounded response."""
    return _prove_bounded_response(
        eps, segCount, data_0, data_1, a, b
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
