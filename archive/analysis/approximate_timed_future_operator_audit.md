# Approximate timed future-operator audit

## Implementation status

The four high-impact items from this audit are now implemented for future
timed `Eventually`, `Always`, and non-strict `Until`:

1. unary truth is computed through four composable Boolean summaries and
   prefix-indexed whole-segment ranges;
2. language concatenation uses a measured hybrid of scalar iteration for
   active spans up to eight and checked bitset shift/OR above that;
3. `possibleUntilBits` uses four RHS-index bitplanes per LHS row with
   word-level transitions and reusable grids; and
4. shifted events are bucketed with monotone sweeps, windows are located by
   index, change bounds use prefix sums, and timed `Until` obtains associative
   whole-range profile products from a segment tree while restricting only
   query boundary segments.

The previous implementations remain callable through `*Legacy` hooks and are
used as exact differential oracles. Eager range-cache construction records
overflow as node state and raises it only when a query selects that range, so
unqueried wide products do not introduce new failures.

Representative end-to-end comparisons against the frozen pre-change header
were:

| Workload | Before | Implemented | Speedup |
| --- | ---: | ---: | ---: |
| `S=128`, rich bounded `[0,8)` Always | 239.7 ms | 1.1 ms | about 220x |
| `S=128`, rich bounded `[0,8)` Until | 664.3 ms | 34.9 ms | 19.0x |
| `S=128`, alternating shifted-unbounded Until | 193.8 ms | 25.2 ms | 7.7x |
| `S=256`, constant shifted-unbounded Until | 324.8 ms | 30.4 ms | 10.7x |

The maintained differential test includes exhaustive small Until-solver
cases, randomized complete operators, closures, clipping, irregular
segmentations, partial ranges, indexed-profile comparisons, and deferred
overflow. The complete CTest suite and ASan/UBSan run pass. The remaining
sections preserve the original baseline and reasoning that motivated the
changes.

## Scope and conclusion

This audit covers only future timed `Eventually`, `Always`, and non-strict
`Until`. Strict future `Until` and strict `Since` are intentionally disabled.
Past timed operators share much of the same machinery, but are not included in
the measurements or the ranking below.

The largest opportunities are not small C++ cleanups:

1. Timed `Always` and `Eventually` construct complete alternating-word
   profiles even though they keep only two truth bits. Replacing those profiles
   with four composable Boolean properties is exact and produced the largest
   measured improvement.
2. Timed `Until` really does need richer profile information. Its hottest
   bounded-input kernel is the scalar Cartesian implementation of
   `bitsetConcat`; an exact word-shift implementation gave a 3.1--3.9x speedup
   on the maintained real inputs.
3. On wide and shifted-unbounded `Until`, the quadratic
   `possibleUntilBits` state space takes over. The word-shift improvement does
   not help that regime; the solver needs its own compact, bit-parallel
   redesign.
4. All three operators repeatedly rediscover shifted endpoints, scan the full
   segmentation for a contiguous window, and rescan ranges to count changes.
   Range indexes and prefix summaries are the next shared structural step.

At the time this audit was first written, the two measured candidate
transformations were kept outside the maintained implementation. They have
since landed behind retained legacy differential hooks, together with the
solver and range-index work described above.

## Pre-change implementation

The shared code is in [`monitoring.hpp`](../../include/monitoring.hpp):

- `bitsetTimedUnary` implements future and past `Always`/`Eventually`.
- `timedProfile` restricts every intersecting segment and concatenates the
  resulting factor languages.
- `bitsetTimedUntil` builds both operand profiles, then calls
  `possibleUntilBits`.
- `bitsetConcat` materializes the language concatenation used for profiles and
  output regions.

For every output segment, the future timed-unary path currently:

1. scans all segmentation endpoints for each shifted window boundary;
2. materializes singleton and open-interval placements;
3. scans all input segments in `timedProfile` for every placement;
4. concatenates the restricted factors into a complete profile;
5. runs a full untimed `Always` or `Eventually` on that profile;
6. keeps only the resulting `can-false` and `can-true` bits;
7. rescans all endpoints and segments through `changes` when an open placement
   is ambiguous;
8. concatenates a structured truth region into the output language.

Timed `Until` performs the same geometry, but constructs two profiles per
placement and feeds them to a reachability search. Its `seen` array has

```text
cells * lhsProfileLength * rhsProfileLength * 4
```

bytes for each active start-value pair. There are at most seven geometric
cells, and up to four start-value pairs. The search is therefore
`Theta(cLR)` in time and space before accounting for profile construction.

The closed, zero-lower-bound unbounded wrappers already route to the untimed
fast paths. The expensive generic path is selected by a positive lower bound,
an open lower endpoint, or a finite upper bound.

## Baseline measurements

All timings below are release-style `-O3` measurements against the current
header. They measure operator execution after input construction. Single-run
numbers are intended to establish order of magnitude and scaling, not stable
benchmark scores.

### Timed Always and Eventually

For unit segments and bounded `[0,8)`, the two unary operators have effectively
the same cost because they share all geometry and profile construction:

| S | Definite Always | Definite Eventually | Ambiguous Always | Ambiguous Eventually |
| ---: | ---: | ---: | ---: | ---: |
| 64 | 2.81 ms | 2.80 ms | 11.67 ms | 10.74 ms |
| 128 | 6.04 ms | 5.93 ms | 23.18 ms | 22.11 ms |
| 256 | 11.05 ms | 12.58 ms | 46.51 ms | 45.84 ms |
| 512 | 25.12 ms | 27.69 ms | 93.85 ms | 93.78 ms |

A shifted-unbounded `[1,infinity)` window exposes the quadratic path even on a
constant language:

| S | Always | Eventually |
| ---: | ---: | ---: |
| 64 | 10.12 ms | 10.11 ms |
| 128 | 39.26 ms | 38.79 ms |
| 256 | 155.86 ms | 160.55 ms |
| 512 | 630.18 ms | 645.42 ms |

Each doubling is approximately 4x.

### Timed Until

On bounded `[0,8)`, timed `Until` costs about 2.2--2.4 times the matching unary
operation for low-complexity languages because it builds two profiles and runs
the solver:

| S | Constant Until | Alternating Until | Ambiguous Until |
| ---: | ---: | ---: | ---: |
| 32 | 2.65 ms | 2.78 ms | 11.16 ms |
| 64 | 5.67 ms | 5.99 ms | 24.71 ms |
| 128 | 11.85 ms | 12.52 ms | 49.91 ms |
| 256 | 24.61 ms | 26.74 ms | 103.28 ms |

Richer factor languages make profile concatenation dominant. At `S = 128`,
bounded `[0,8)` took 248.7 ms for timed `Always` and 676.8 ms for timed
`Until`. The `Until` call requested about 75.7 MB across 18,000 allocations.

An instrumented `S = 64` rich-language run attributed 90.3% of sampled time to
`bitsetConcat` and 9.7% to `possibleUntilBits`. This profile combines the two
unary calls and one `Until` call, but clearly identifies the bounded rich-input
hot path.

Shifted-unbounded alternating inputs move the bottleneck into the solver:

| S | Timed Until | Requested bytes/call |
| ---: | ---: | ---: |
| 64 | 34.99 ms | 7.57 MB |
| 128 | 193.02 ms | 40.97 MB |
| 256 | 1,219.58 ms | 252.39 MB |

The growth is worse than the already-quadratic unary geometry because profile
lengths grow toward the trace end and the solver is quadratic in the two
profile maxima.

### Isolated Until solver

With one length-`L` word in each operand and an open window, direct
`possibleUntilBits` timings were:

| L | Time range over start-value pairs | Requested bytes |
| ---: | ---: | ---: |
| 16 | 16.5--17.6 us | 4.2 KB |
| 32 | 44.7--57.4 us | 14.4 KB |
| 64 | 217--222 us | 53.3 KB |
| 128 | 681--922 us | 204.9 KB |
| 256 | 2.80--3.72 ms | 802.9 KB |
| 512 | 11.43--11.54 ms | 3.18 MB |

Both time and the `seen` allocation follow the expected quadratic curve.

## Validated candidate 1: Boolean timed-unary summaries

For any nonempty profile language `P`, timed unary needs only:

| Property | Exact profile predicate |
| --- | --- |
| `allZero` | `P[0][0]` |
| `allOne` | `P[1][0]` |
| `hasZero` | `P[0].any()` or a non-singleton word in `P[1]` |
| `hasOne` | `P[1].any()` or a non-singleton word in `P[0]` |

These properties compose exactly over factor concatenation:

```text
allZero(L ++ R) = allZero(L) AND allZero(R)
allOne (L ++ R) = allOne (L) AND allOne (R)
hasZero(L ++ R) = hasZero(L) OR  hasZero(R)
hasOne (L ++ R) = hasOne (L) OR  hasOne (R)
```

The final unary truth bits are:

```text
Always:     canFalse = hasZero, canTrue = allOne
Eventually: canFalse = allZero, canTrue = hasOne
```

A temporary full-operator candidate retained the current geometry,
restriction helpers, `changes`, and output composition, replacing only the
materialized profile and untimed reduction. It matched the current output in
20,000 randomized cases covering:

- bounded and shifted-unbounded windows;
- all lower/upper closure combinations;
- punctual, empty-after-clipping, and ordinary windows;
- irregular segmentations;
- both unary operators; and
- partial `[s,e)` evaluation.

At `S = 128`, bounded `[0,8)` kernel timings were:

| Language | Current Always | Summary Always | Speedup |
| --- | ---: | ---: | ---: |
| Constant | 5.51 ms | 1.01 ms | 5.5x |
| Ambiguous | 24.09 ms | 1.44 ms | 16.7x |
| Rich | 248.62 ms | 1.47 ms | 169.6x |

`Eventually` produced the same pattern.

On the maintained real signal inputs, timed `Always` produced:

| Samples / S | Window | Current | Summary | Speedup |
| ---: | ---: | ---: | ---: | ---: |
| 64 / 60 | 1 s | 3.56 ms | 2.16 ms | 1.64x |
| 64 / 60 | 2 s | 6.80 ms | 1.96 ms | 3.46x |
| 128 / 127 | 1 s | 6.82 ms | 3.31 ms | 2.06x |
| 128 / 127 | 2 s | 14.37 ms | 3.32 ms | 4.33x |
| 256 / 249 | 1 s | 15.04 ms | 7.40 ms | 2.03x |
| 256 / 249 | 2 s | 32.47 ms | 6.91 ms | 4.70x |
| 512 / 477 | 1 s | 23.73 ms | 13.09 ms | 1.81x |
| 512 / 477 | 2 s | 52.81 ms | 13.14 ms | 4.02x |

The summary also prevents overflow caused solely by constructing an
intermediate unary profile. The existing final output-region capacity check
must remain.

The temporary candidate intentionally retained full segmentation scans and
allocated restricted factors. On constant shifted-unbounded `S = 512`, it
still requested 67.8 MB, despite reducing runtime from 636.8 ms to 8.34 ms.
The production version should therefore combine the Boolean algebra with the
range-index work below rather than stop at the prototype.

## Validated candidate 2: word-shift concatenation for Until

The current scalar `bitsetConcat` visits every right index for every active
left index. For a fixed left endpoint `i`, its entire contribution from a
right bucket is exactly that bucket shifted by

```text
i + (lastLeftValue != rightFirstValue).
```

The candidate replaces the scalar right-index loop with a checked
`std::bitset` shift and OR. It preserves the existing capacity failure by
checking the right maximum before shifting.

Validation covered:

- 20,000 randomized concatenation pairs; and
- 10,000 complete timed-`Until` evaluations over bounded/unbounded windows,
  all closures, irregular segmentations, and partial ranges.

Every result and overflow outcome matched the current implementation.

On the maintained real signal inputs:

| Samples / S | Window | Current Until | Shift Until | Speedup |
| ---: | ---: | ---: | ---: | ---: |
| 64 / 60 | 1 s | 8.99 ms | 2.92 ms | 3.08x |
| 64 / 60 | 2 s | 17.35 ms | 4.56 ms | 3.81x |
| 128 / 127 | 1 s | 22.06 ms | 6.85 ms | 3.22x |
| 128 / 127 | 2 s | 46.03 ms | 11.93 ms | 3.86x |
| 256 / 249 | 1 s | 47.05 ms | 14.88 ms | 3.16x |
| 256 / 249 | 2 s | 98.22 ms | 25.21 ms | 3.90x |

At `S = 128`, bounded `[0,8)` improved from 51.0 to 13.9 ms on an
ambiguous language and from 671.4 to 143.3 ms on a rich language. Sparse
constant and single-word alternating inputs were unchanged because shifting a
full fixed-capacity bitset is not cheaper than scanning a tiny right language.
A production kernel should retain the scalar path for very sparse right
buckets and select the shift path by density or active-index count.

The word-shift candidate did not improve shifted-unbounded alternating Until:

```text
S=64:   34.87 ms current, 34.97 ms shift
S=128: 194.22 ms current, 192.76 ms shift
S=256:   1.208 s current,  1.207 s shift
```

That is the regime where `possibleUntilBits`, not concatenation, dominates.

## Original ranked implementation plan

### 1. Replace timed-unary profiles with range-composable truth summaries

Potential impact: **highest and already demonstrated**.

Implement the validated four-Boolean algebra for future `Always` and
`Eventually`. Precompute whole-segment summaries and prefix counts, and compute
only the one or two restricted boundary summaries per window. Keep the old
path as a test oracle until the full differential matrix passes.

This removes width-dependent profile growth, most hot allocations, calls into
untimed unary, and intermediate-profile capacity failures. Measured gains are
1.6--4.7x on current real workloads and much larger on complex languages.

### 2. Add a hybrid word-shift `bitsetConcat` path for timed Until

Potential impact: **high on current real bounded Until**.

Use the validated shift/OR formula when the right factor is sufficiently
dense; retain the current sparse scalar loop for singleton and low-popcount
factors. Hoist both right maxima, keep checked overflow, and apply the hybrid
kernel to profile and output-region concatenation.

The all-shift prototype already gives 3.1--3.9x on maintained real inputs.
The hybrid selection avoids its neutral or slightly negative sparse cases.

### 3. Redesign `possibleUntilBits` as a compact bit-parallel reachability solver

Potential impact: **highest for wide or shifted-unbounded Until**.

As an incremental first step, replace four `char` entries per
`(cell,lhsIndex,rhsIndex)` with one four-bit state mask and reuse buffers across
placements. This cuts the dominant `seen` storage by 4x and removes repeated
large allocations.

The main rewrite should process RHS-index reachability as machine-word
bitsets for each LHS row and flag state, exploiting the monotone grid. The
target bound is approximately `O(cLR/word)` time with row- or cell-local
storage, while preserving exact terminal membership checks against the two
profile bitsets.

This work is independent of word-shift concatenation: the latter produced no
speedup at `S = 256` shifted-unbounded Until, where the current call takes
about 1.21 s and requests 252 MB.

### 4. Index contiguous windows and shifted critical events

Potential impact: **high asymptotically and shared by all three operators**.

Apply three related changes:

1. locate the first and last intersecting segments with binary search or
   monotone pointers instead of scanning all `S` segments in every profile or
   summary query;
2. bucket each shifted endpoint into its containing output segment once,
   replacing `Theta(kS^2)` critical-candidate checks with `O(kS)` events plus
   placement output; and
3. for timed Until, cache associative whole-range profile products in a range
   tree or sliding-window monoid, leaving only restricted boundary factors to
   concatenate per placement.

Unary summary range queries can become constant time for full interiors.
Until range products remain capacity-bounded but avoid reconstructing the same
overlapping products for every moving window.

### 5. Prefix-sum the `changes` bound

Potential impact: **medium on ambiguity-heavy inputs**.

Precompute each segment's
`max(msb(bucket0), msb(bucket1))` weight and its saturated prefix sum. Count
internal endpoints with two bounds and answer the intersecting segment-weight
sum as a range query. Timed Until needs one prefix structure per operand.

This replaces every ambiguous placement's endpoint and segment rescans with
`O(log S)` searches or monotone-pointer `O(1)` queries. Saturation at `N`
preserves the current output-capacity exception without risking integer
overflow.

### 6. Specialize output-region composition

Potential impact: **medium after the larger profile work**.

Every placement emits only one of:

- constant false;
- constant true;
- both singleton values; or
- equal dense prefixes in both buckets.

Maintain the output in an endpoint/extrema representation or compose these
four region forms with specialized shift/OR operations. Avoid sending them
through the fully generic Cartesian concatenator.

### 7. Remove wrapper copies and fixed-small allocations

Potential impact: **low to medium; useful after algorithmic work**.

Change the bounded future wrappers to take languages and segmentation by const
reference. Use `array` for at most three offsets, reserve critical and
placement buffers, reuse scratch storage across output segments, and avoid
allocating two-bucket vectors in inner loops.

These changes matter on sparse short-window inputs but do not alter the
quadratic geometry or solver.

### 8. Add validation and production differential hooks

Potential impact: **correctness/risk reduction rather than speed**.

Before switching implementations, validate equal operand sizes, exactly
`S + 1` strictly increasing segmentation endpoints, nonnegative ordered
bounds, two buckets per segment, and valid `[s,e)` ranges. Keep callable old
implementations in tests and compare complete bitsets, not only verdicts.

The matrix should include all four closure pairs, punctual windows, windows
outside the domain, clipping at both domain endpoints, irregular
segmentations, partial ranges, empty/near-capacity results, and overflow
equivalence. Timed Until also needs isolated solver tests at high profile
lengths and every geometric cell signature.

## Recommended implementation sequence

The practical sequence is:

```text
unary differential hook
  -> four-Boolean timed Always/Eventually
  -> indexed boundary/range summaries
  -> hybrid shift/OR concatenation
  -> timed-Until differential hook
  -> compact then bit-parallel Until solver
  -> shifted-event buckets and changes prefix sums
  -> specialized output composition and API cleanup
```

This sequence lands the strongest measured unary improvement first, then the
strongest measured bounded-Until improvement, while establishing the oracle
coverage needed for the more invasive wide-window solver rewrite.

## Past-time port implementation results

The optimized timed infrastructure now also serves bounded and shifted-
unbounded past `Once`, `Historically`, and non-strict `Since`. The port keeps
the future paths intact and introduces an explicit `TimedDirection` instead
of duplicating subtly different geometry code.

### Implemented architecture

- Direction-aware monotone critical-event bucketing uses `endpoint - offset`
  for future windows and `endpoint + offset` for past windows. The shared
  change index queries `(left + offset, right + offset)` in the future and
  `(left - offset, right - offset)` in the past.
- Timed `Once` and `Historically` use the same four composable truth summaries
  as their future counterparts: `allZero`, `allOne`, `hasZero`, and `hasOne`.
  Their truth reductions are direction-independent, so past monitoring no
  longer materializes an alternating profile merely to discard all but two
  truth possibilities.
  `Historically` uses `canFalse = hasZero` and `canTrue = allOne`;
  `Once` uses `canFalse = allZero` and `canTrue = hasOne`.
- Timed non-strict `Since` builds one associative `TimedProfileIndex` and one
  prefix `TimedChangeIndex` per operand. It queries each horizon in
  chronological order, reverses the completed operand profiles and timed
  ranges, and reuses the packed `possibleUntilBits` solver. Output regions are
  assembled chronologically with the hybrid scalar/word-shift concatenator.
- Empty factors retain a complete-call fallback because cached associative
  profile products require nonempty factors. Cached overflow is deferred and
  raised only when a query selects the overflowing node.
- All shifted-coordinate arithmetic is widened before adding, subtracting,
  doubling, or negating. Shared validation rejects malformed shapes,
  segmentation, intervals, and ranges deterministically.

`bitsetTimedUnaryLegacy` and the fully independent
`bitsetTimedSinceLegacy` remain differential hooks. The Since hook is pinned
to legacy profile construction, scalar Until reachability, scalar
concatenation, and nested geometry scans, so it cannot accidentally inherit
the optimized components it is intended to check.

### Exact past-time geometry

For doubled evaluation coordinate `t`, the implementation preserves:

| Operation | Finite range | Upper-infinite range | Closures |
| --- | --- | --- | --- |
| Once/Historically window | `[t-2b, t-2a]` | `[domain.left, t-2a]` | finite `{rightClosed,leftClosed}`; infinite `{true,leftClosed}` |
| Since horizon | `[t-2b, t]` | `[domain.left, t]` | closed/closed |
| Since witness window | `[t-2b, t-2a]` | `[domain.left, t-2a]` | finite `{rightClosed,leftClosed}`; infinite `{true,leftClosed}` |

All ranges are intersected with the existing half-open timed domain
`[0, 2 * segmentation.back())`. At `a == b`, only closed/closed is a
nonempty punctual interval. Empty-window vacuity, the punctual
`segmentFirstBit(timedProfile(...))` rule, witness-before-LHS ordering, output
capacity checks, and zero padding outside `[s,e)` are unchanged.

### Focused verification

The new timed-past target compares complete languages and exception classes.
It includes:

- 20,000 randomized direction-specific critical-bucket cases, each with five
  past change-index queries per distinct offset, plus near-limit arithmetic;
- 6,000 randomized timed-unary inputs, each checked for both `Once` and
  `Historically`, giving 12,000 primary optimized/legacy operator comparisons;
- 12,000 randomized timed-Since comparisons;
- all bounded and upper-infinite interval families, all closure pairs,
  punctual and clipped windows, partial ranges, retained history, empty-factor
  fallback, and input immutability; and
- fixed endpoint/closure cases, range-reversal involution, word boundaries
  63/64/127/128, wrapper routing, and malformed calls.

The complete Release CTest run passed 7/7. The focused ASan/UBSan run passed
4/4 selected regression executables, including the existing timed-future
suite. Strict Since remains disabled under `#if 0` and is still absent from
the maintained target and dispatcher lists.

### Measured performance

All values below are Release/O3 median speedups from seven samples after one
warm-up, with full-language equality checked before timing. `Ambiguous` means
each factor contains the false singleton and true singleton. `Rich` means both
singletons plus one length-two and one length-three alternating word; LHS and
RHS use opposite phases.

For bounded `[0,8)` at 256 and 512 segments, speedup over the fully frozen
legacy implementation was:

| Operator | Constant | Ambiguous | Rich |
| --- | ---: | ---: | ---: |
| Timed Once/Historically | 5.94--6.14x | 23.3--23.7x | 46.1--46.9x |
| Timed non-strict Since | 1.56--1.62x | 5.74--5.83x | 10.82--11.08x |

For timed Since alone, the final indexed implementation was also 1.41--1.59x
faster at 256/512 segments across all three patterns than the frozen
pre-index implementation that already used the packed solver and hybrid
concatenation. This isolates the benefit of the past-direction event, change,
and profile indexes from the earlier shared optimizations.

Shifted `[1,infinity)` exposes the avoided repeated profile construction more
strongly:

| Fixture | Timed Once/Historically | Timed Since |
| --- | ---: | ---: |
| 256 segments, constant | about 120x | 13.83x |
| 64 segments, ambiguous | about 625--630x | 36.18x |
| 64 segments, rich | about 1,803--1,805x | 42.91x |

The 256-segment ambiguous shifted-Since legacy run exceeded two minutes and
was terminated, so the reported ambiguous and rich shifted comparisons use
64 segments. These are controlled synthetic operator measurements; no
real-workload speedup is claimed here.
