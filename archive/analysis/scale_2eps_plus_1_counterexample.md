# Counterexample to `scale = 2 * eps + 1`

## Verdict

`scale = 2 * eps + 1` is **not sufficient** for exact timed-response
verdicts. The failure occurs in both the ideal inverse-knot grid model and the
literal cell-based encoding in
`benchmarks/offline/exact/smt/random.py`.

The counterexample is a completeness failure for the satisfying polarity: a
semantically admissible continuous retiming satisfies the formula, but the
scale-3 grid contains no satisfying assignment.

## Instance

Use

```text
d = 9
eps = 1
[a,b) = [3,4)
p = 100100000    terminal value: 0
q = 000010001    terminal value: 1
```

The formula is

```text
G(p -> F_[3,4) q).
```

Since `eps = 1`, the proposed scale is `K = 2*eps + 1 = 3`.

## Exact continuous witness

Let `R_j = j` be the reference inverse clock. Define the two agent inverse
clocks at the integer local boundaries by

| `j` | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `P_j` | 0 | 1/4 | 5/4 | 7/2 | 15/4 | 17/4 | 21/4 | 25/4 | 29/4 | 9 |
| `Q_j` | 0 | 1/4 | 5/4 | 9/4 | 15/4 | 17/4 | 21/4 | 25/4 | 29/4 | 9 |

Interpolate linearly between consecutive knots. Each row is strictly
increasing. Moreover, for every `j = 0,...,8`,

```text
max(R_j, P_j, Q_j) < min(R_(j+1), P_(j+1), Q_(j+1)).
```

This implies every ordered pair of inverse clocks satisfies
`A(x) < B(x+1)` for all real `x` where it is defined: on each unit interval,
the difference is linear and is strictly negative at both endpoints.
Consequently, the reconstructed forward clocks have strict pairwise skew
less than `eps = 1` at every physical time. The same inequalities preserve
the mandatory happened-before order.

Under this retiming, the true physical-time intervals are

```text
p: [0, 1/4) union [7/2, 15/4)
q: [15/4, 17/4) union [29/4, 9).
```

A `q` interval `[S,T)` responds to exactly the trigger times in the open
window `(S-b,T-a)`. Here the two response windows are

```text
(-1/4, 5/4) and (13/4, 6).
```

They respectively contain `[0,1/4)` and `[7/2,15/4)`. Therefore the dense
retiming satisfies `G(p -> F_[3,4) q)`.

The identity retiming violates the formula: at trigger time `u = 0`, `q` is
false throughout `[3,4)`. Thus the semantic verdict is inconclusive, not
definitely false.

## Why the scale-3 knot grid cannot satisfy the formula

Write `P_j` and `Q_j` for arbitrary admissible inverse knots. Any satisfying
retiming of this instance must obey

```text
3 < Q_8 - 4 < P_3 < Q_4 < 4.                 (1)
```

The four parts are forced as follows.

1. Strict skew against the reference gives `7 = R_7 < Q_8`, hence
   `3 < Q_8 - 4`.
2. The first `q` response window ends at `Q_5 - 3 < R_6 - 3 = 3`, while the
   second starts at `Q_8 - 4 > 3`. The windows are separated around time 3.
3. Strict skew gives `R_3 < P_4`, so the second `p` interval `[P_3,P_4)`
   extends beyond time 3. It cannot be covered by the first response window
   or cross the gap between the windows. It must therefore start inside the
   second window, giving `Q_8 - 4 < P_3`.
4. Strict cross-clock skew gives `P_3 < Q_4`. At trigger `u = 0`, only the
   first `q` run can respond, so its start must satisfy `Q_4 < 4`.

At scale `K = 3`, all three middle quantities in (1) are multiples of `1/3`.
There are only two grid points strictly between 3 and 4, namely `10/3` and
`11/3`. Three strictly ordered grid values cannot fit there. Hence no
scale-3 satisfying knot model exists.

## Why the literal cell encoding also fails

For the literal encoding, let `P_j` and `Q_j` denote the first integer ticks
at which the corresponding monotone, surjective cell clock enters local cell
`j`. The strict padded-domain predicates give

```text
(j-1)K < P_j,Q_j < (j+1)K
```

for `eps = 1`. The happened-before predicate gives `P_3 < Q_4`. Applying the
same connected-window argument directly to the quantified tick query forces

```text
3K < Q_8 - 4K < P_3 < Q_4 < 4K.              (2)
```

For `K = 3`, (2) requires three distinct integer ticks strictly between 9
and 12, but only 10 and 11 are available. Therefore the literal maintained
encoding also has no satisfying assignment.

## Computational check

The checks used Z3 4.16.0.

```text
archived optimized continuous solver:
    positive=False, negative=False
    verdict=inconclusive

literal maintained encoding with scale patched in memory to K=3:
    positive=False, negative=True
    verdict=false

literal maintained encoding with K=4:
    positive=False, negative=False
    verdict=inconclusive
```

An independent ideal inverse-knot encoding likewise returned

```text
continuous rationals: sat
K=3 integers:         unsat
K=4 integers:         sat
```

The scale-4 satisfying knot assignment is exactly four times the rational
witness above.

## Consequence

Do not replace the current scale by `2*eps + 1`: it still produces an
incorrect definite-false verdict on this instance. `K = 4` repairs this
particular counterexample but is not a universal bound. The currently proved
general grid fallback is the duration-dependent scale `K = 2*d - 1`; finding
a smaller universally sufficient bound remains a separate problem.

