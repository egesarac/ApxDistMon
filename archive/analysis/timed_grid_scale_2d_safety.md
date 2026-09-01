# Why timed grid scales `2*d - 1` and `2*d` are safe

## Result

For the two-agent, duration-`d` timed-response problem, the inverse-clock
difference-constraint model is complete on a grid with

```text
K = 2*d - 1
```

ticks per physical-time unit. Consequently, the slightly larger and simpler
choice

```text
K = 2*d
```

is also safe.

Here, *safe* means that discretization preserves the existence of both:

1. an admissible retiming satisfying the timed formula; and
2. an admissible retiming violating the timed formula.

Preserving both existential polarities preserves the three-valued verdict.
These are sufficient upper bounds, not claims of minimality.

## Scope and assumptions

The argument assumes:

- two agents;
- integer duration `d >= 1`;
- integer skew and temporal bounds;
- Boolean signals whose changes occur at integer local times;
- strictly increasing inverse clocks;
- strict skew and happened-before constraints;
- a timed-response formula whose dense-time truth and falsity have been
  reduced to finite branches over signal-run endpoints.

The maintained timed formula

```text
G(p -> F_[a,b) q)
```

has such a reduction. For a satisfying branch, each true `p` interval must
be covered by a connected chain of shifted true-`q` windows. For a violating
branch, some true `p` interval is not covered. After choosing the relevant
windows or uncovered gap, no additional real-time witness variable is
necessary.

## Counting the variables

For each agent, introduce inverse-clock boundary values

```text
z[i,0], z[i,1], ..., z[i,d].
```

The endpoint values are fixed:

```text
z[i,0] = 0
z[i,d] = d.
```

Therefore each agent has only `d-1` free knots. With two agents, there are

```text
2*(d-1)
```

free variables. Add one constant-zero anchor used to represent integer
constants in the constraint graph. The graph therefore has

```text
N = 2*(d-1) + 1 = 2*d - 1
```

vertices.

The reference clock is fixed and does not introduce free variables.
Expressions such as `Q_j-a` and `Q_j-b` are translations of the same knot,
not new variables.

## Finite formula branches are strict difference constraints

Once a Boolean formula branch is fixed, every condition has one of the forms

```text
x_v - x_u <= c
x_v - x_u <  c
```

with integer `c`. This covers:

- inverse-clock monotonicity;
- strict shifted inequalities implementing clock skew;
- happened-before constraints;
- response-window containment;
- overlap between consecutive response windows;
- finite-horizon and endpoint conditions.

An equality is represented by two non-strict inequalities. Thus a branch is
a strict difference-bound system, or strict DBM.

## The tighter bound: `K = 2*d - 1`

Consider a directed cycle `C` in the branch's difference-constraint graph.
Let:

- `W(C)` be the sum of its integer edge constants;
- `s(C)` be the number of strict edges in the cycle.

Dense feasibility implies

```text
W(C) >= 0.
```

If the cycle contains a strict edge, dense feasibility requires

```text
W(C) > 0.
```

Otherwise summing the cycle would yield `x < x`. Since `W(C)` is an
integer,

```text
s(C) > 0  implies  W(C) >= 1.                 (1)
```

Now use a scale-`K` grid and write `X = K*x`, where `X` is integral.
A non-strict constraint becomes

```text
X_v - X_u <= K*c,
```

while a strict constraint becomes

```text
X_v - X_u <= K*c - 1.
```

The resulting integer weight around cycle `C` is therefore

```text
K*W(C) - s(C).                                (2)
```

It is enough to consider simple cycles: every negative closed walk contains
a negative simple cycle. A simple cycle visits each graph vertex at most once,
so

```text
s(C) <= length(C) <= N = 2*d - 1.             (3)
```

Choose `K=N=2*d-1`. If `C` has no strict edge, (2) is nonnegative because
`W(C)>=0`. If it has a strict edge, combine (1)--(3):

```text
K*W(C) - s(C)
    >= N*1 - N
    = 0.
```

Thus scaling creates no negative cycle. An integer-weight difference system
without a negative cycle has an integer solution, obtained for example from
shortest-path distances. Dividing that solution by `K` gives a grid retiming
that satisfies the same branch.

Therefore every dense-feasible formula branch has a representative on the
`K=2*d-1` grid.

The reverse direction is immediate for the knot model: a grid assignment is
a rational dense assignment satisfying the same inequalities. Hence the grid
is sound as well as complete.

## The simpler bound: `K = 2*d`

The scale

```text
K = 2*d
```

is safe by the same cycle argument because

```text
2*d >= 2*d - 1.
```

This is not an appeal to grid nesting: denominator-`2*d-1` points need not
belong to the denominator-`2*d` grid. Rather, equations (1)--(3) prove
directly that **every** integer `K >= N` is sufficient, because
`K*W(C)-s(C) >= K-N >= 0` on every strict simple cycle.

It also has a simple slot interpretation. A formula branch has only
`N=2*d-1` independent graph levels, counting the constant anchor. A simple
strict dependency chain cannot use more than these `N` levels before
repeating a vertex.

A scale-`2*d` grid has exactly `2*d-1` grid points strictly inside each
integer interval:

```text
1/(2*d), 2/(2*d), ..., (2*d-1)/(2*d).
```

Even in the worst case, there is therefore one interior slot for every
strictly ordered level that the branch can force across one unit of integral
slack. Non-strictly related levels may share a slot.

This slot explanation is the intuitive version of the cycle proof. The cycle
proof remains the formal justification: replacing `K=2*d-1` by the larger
`K=2*d` only increases every positive term `K*W(C)`, so it cannot create a
negative cycle.

## Why the proof covers the full verdict

The timed formula is disjunctive: different retimings can select different
response-window chains or different violating gaps. The proof does not require
one grid assignment to preserve every branch simultaneously.

Instead:

1. choose any dense-feasible satisfying or violating branch;
2. apply the difference-constraint argument to that branch;
3. obtain a grid assignment satisfying the same branch.

Therefore:

```text
some dense satisfying retiming exists
iff
some grid satisfying retiming exists,
```

and independently:

```text
some dense violating retiming exists
iff
some grid violating retiming exists.
```

This is exactly what is required for the positive, negative, and inconclusive
verdicts.

## Relation to the maintained cell-based encoding

The proof is stated using inverse-clock entry knots. An integer knot assignment
induces a cell clock by defining

```text
c_i(t) = max { j : z[i,j] <= t }.
```

Thus `z[i,j]` is the first grid tick at which clock `i` enters local cell
`j`. The resulting clock is monotone and surjective. The shifted knot
inequalities imply the maintained padded cell-domain constraints, and the
variation-knot ordering implies the maintained happened-before constraints.
The dense response branch becomes the corresponding quantified tick branch.

Conversely, first-entry ticks extracted from a literal cell assignment provide
the relevant inverse knots. At signal variations, happened-before supplies the
required strict ordering. Equalities at silent boundaries can be strictified
without changing either signal value or formula polarity.

Accordingly, the knot completeness bound supplies the intended safe scale for
the maintained grid representation as well.

## Cost of using `2*d` instead of `2*d - 1`

The maintained implementation explicitly generates a quadratic family of
monotonicity constraints over a horizon of `K*d` ticks. Its dominant
constraint-size proxy is therefore

```text
O((K*d)^2).
```

The estimated overhead of `2*d` relative to `2*d-1` is

```text
(2*d / (2*d-1))^2.
```

Examples:

| `d` | estimated `2*d` versus `2*d-1` |
|---:|---:|
| 4 | 1.31x |
| 8 | 1.14x |
| 16 | 1.07x |

The relative overhead tends to one as `d` grows. Thus `2*d` is a reasonable
choice when a simpler statement and implementation are preferred over the
small tightening obtained from `2*d-1`.

## Limitations

- Neither bound is generally minimal.
- The counterexample ladder shows that no bound depending only on `eps` can
  be universally sufficient.
- For `eps=1`, that family requires scales growing approximately as
  `(d-1)/2`, leaving a substantial gap below the upper bounds here.
- The vertex count is specific to two agents. With more agents, the same proof
  must use the corresponding larger number of free inverse knots.
