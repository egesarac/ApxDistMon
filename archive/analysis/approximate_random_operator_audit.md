# Approximate Monitoring Efficiency Audit: AND/OR and Always/Eventually

Date: 2026-08-22

## Scope

This report is intentionally limited to these approximate-monitor operations:

- conjunction: `bitsetConjunction`;
- disjunction: `bitsetDisjunction`;
- untimed future always/eventually: `bitsetAlways` and
  `bitsetEventually`;
- bounded and unbounded future always/eventually through
  `bitsetTimedUnary`.

The focus is the random-signal workload. Other temporal operators, formula
composition, preprocessing performance, case-study workloads, exact
monitoring, and the online algorithm architecture are outside this report.

The checked-in offline random driver normally stops at 32 samples
([random.cpp:255](../../benchmarks/offline/approximate/random.cpp#L255)). To expose
scaling, I also used 64, 128, 256, and 512 samples, which are normal online
durations ([run_online.py:71](../../scripts/run_online.py#L71)), a real
1,024-sample stress input, and controlled inputs through 1,024 segments. Signal
conversion, uncertainty construction, canonical segmentation,
value-expression construction, and AP conversion were all completed before
the operator timers started.

No monitoring implementation was changed. All prototypes, allocation probes,
and microbenchmarks were built under `/tmp`; this file is the only repository
artifact produced by this audit.

The current random driver reaches disjunction and eventuality through explicit
Boolean duality rather than calling the standalone `bitsetDisjunction` or
`bitsetEventually` functions. Standalone-helper measurements below are
operator/library evidence. The expanded OR shell is separately identified,
and the eventually results affect the current workload only if evaluation
adopts a direct eventual path.

## Bottom line

There are real and substantial inefficiencies, but they differ by operator:

| Operator path | Main cause | Long-input observation | Best direction |
| --- | --- | ---: | --- |
| Conjunction | Eight scalar reverse searches through 1,000-bit sets; by-value deep copies | About 4.1--4.2 us/segment | One word-level extrema pass, const inputs, flat segments |
| Disjunction | Same conjunction cost plus materialized De Morgan negations and copies | Expanded random-driver shell about 4.1--4.2 us/segment; standalone lvalue helper about 4.2--4.3 | Polarity-aware direct binary kernel |
| Untimed always | Computes four maxima although it needs only three existence predicates | 2.092 ms at 1,024 active segments; direct predicates 25.0 us | Replace maxima with exact predicates |
| Untimed eventually | Same wasted maxima, with the dual carry condition | 2.094 ms at 1,024 active segments; direct predicates 27.4 us | Replace maxima with exact predicates |
| Timed always/eventually | Repeated all-segment searches, full profile construction, generic concatenation, repeated range scans | Actual 2-second bound about 55 ms at 477 segments; shifted-unbounded and wider stress paths reach about 648/618 ms | Evaluate four Boolean range summaries instead of constructing profiles |

The most immediately attractive change is the validated untimed predicate
rewrite: approximately 70--84x faster when the carry remains active on the
low-index languages characteristic of these random inputs. The largest
algorithmic problem is timed unary evaluation: its geometry performs
quadratic full-segmentation scans, and its profile construction adds a
window-width-dependent explosion on top.

There is also a correctness defect in conjunction: it can compute an index
greater than 999 and write it with unchecked `bitset::operator[]`.

## Measurement method and notation

Measurements used GCC 13.3, C++17, release-style `-O3 -DNDEBUG`, on an Intel
Core Ultra 7 255U. The untimed microbenchmarks additionally used
`-march=native`. Each comparison in a table was compiled into the same
binary; absolute times should not be compared across different tables.
Checksums or sampled output bits prevented dead-code elimination.

Random languages came from `data/signals/<duration>_0.txt` and
`<duration>_100.txt`, with an 8-second uncertainty bound. The first
64--512-sample pair produced:

| Samples | Canonical segments | AP 0 max index | AP 1 max index | AND max index | OR max index |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 64 | 60 | 12 | 7 | 18 | 18 |
| 128 | 127 | 10 | 12 | 20 | 20 |
| 256 | 249 | 13 | 12 | 24 | 24 |
| 512 | 477 | 13 | 13 | 24 | 25 |

Thus the AP inputs occupy at most the bottom 1.4% of the fixed bitset, and
these AND/OR intermediates at most the bottom 2.6%.
Controlled binary tests take the first `S` random-derived segments below
477, then cycle through them when `S > 477`, reaching exact sizes
64--1,024. Cycling does not change these pointwise binary operators.

Notation used below:

- `S`: number of canonical segments;
- `N`: bitset capacity, currently 1,000;
- `A`: number of segments on which an untimed carry keeps the expensive
  branch active;
- `k`: number of distinct timed-window offsets, at most two here;
- `C`: number of timed placements/windows evaluated;
- `w`: bounded-window width in unit-segment experiments.

## Representation costs shared by all four operators

The capacity is fixed by
[`SIZE = 1000`](../../include/monitoring.hpp#L28). A segment language is a
two-element vector of `bitset<1000>`; a complete language is effectively:

```cpp
vector<vector<bitset<1000>>> language;
```

Bucket 0 contains alternating words beginning with 0, bucket 1 those beginning
with 1, and bit index `i` represents a word of length `i + 1`. Consequently:

| Bucket | Index parity | Endpoint class |
| ---: | ---: | --- |
| 0 | even | `0...0` |
| 0 | odd | `0...1` |
| 1 | odd | `1...0` |
| 1 | even | `1...1` |

On this libstdc++ ABI, `sizeof(bitset<1000>) == 128`. One nested language
therefore requests approximately `280S` bytes:

- `24S` bytes for the outer vector's inner-vector objects;
- `256S` bytes for two bitsets per inner allocation;
- `S + 1` allocations: one outer block and one inner block per segment.

The nested vector is a persistent source of allocation calls, pointer chasing,
and poor locality. A `vector<array<bitset<N>, 2>>` needs one allocation and
`256S` bytes.

The generic [`msb`](../../include/monitoring.hpp#L51) is another shared problem.
It reads one proxy bit at a time from `N - 1` downward. It does not inspect
64-bit words or use a count-leading-zeros instruction. With maximum index 13,
one call tests about 986 irrelevant high zero positions before reaching useful
data. The current expressions combine global mutable parity masks at
[monitoring.hpp:30](../../include/monitoring.hpp#L30) with each input, forming
full masked bitset temporaries before this scalar search.

The functions are templated on `N`, but their masks have the unrelated fixed
type `bitset<SIZE>`, so they are not genuinely capacity-generic. Merely
shrinking `SIZE` would make the random inputs faster, but is not a sound fix:
binary and timed composition can require much longer words, and the timed
stress run already exceeds 1,000. The useful representation needs an active
maximum or word-level metadata rather than paying for a worst-case global
capacity on every search.

### Capacity required at `d = 128`

`SIZE` is a capacity: an active maximum index of `k` requires `SIZE >= k + 1`.
For the 128-sample, 1 Hz random workload with `epsilon <= 8`, the tight and
observed requirements are:

| scope | largest index | minimum `SIZE` |
|---|---:|---:|
| fixed 100-pair corpus, untimed operators | 25 | 26 |
| all generator-valid inputs, untimed operators | 32 | 33 |
| fixed 100-pair corpus, including timed `phi5`/`phi6` | 219 | 220 |
| all generator-valid inputs, timed `phi5` | 278 | 279 |
| all generator-valid inputs, timed `phi6` | 339 | 340 |

For the untimed bound, conversion can place at most index 16 in either AP:
at 1 Hz, at most 16 edge-uncertainty regions overlap one canonical segment.
Conjunction adds operand indices, so both direct AND and the conjunction used
by De Morgan OR can reach `16 + 16 = 32`. Opposite-phase inputs that alternate
at every sample attain this bound.

For timed `phi5` (`[0,1s)`) and `phi6` (`[0,2s)`), the shifted-endpoint
geometry permits at most 9 and 10 critical points. The tight correlated
`changes` sums are 126 and 155; output concatenation therefore reaches inner
indices `2*126 + 17 = 269` and `2*155 + 19 = 329`. The following conjunction
with `p` adds at most 9 or 10, giving indices 278 and 339. Generator-valid
witnesses are constant until samples 119 and 118 respectively, then alternate
through sample 127; the transient profile maxima (152 and 169) are smaller.
The fixed-corpus figures are the maxima from scanning all 100 configured
signal pairs and all four epsilon values.

These capacities are deliberately workload-specific. They cover the random
benchmark's AND/OR, untimed always/eventually, and timed `phi5`/`phi6` paths;
they do not cover `phi3`/until, other durations or sampling rates, or any case
study input.

## Conjunction

### What the current kernel does

[`bitsetConjunction`](../../include/monitoring.hpp#L4164) takes both complete
languages by value because it clears singleton bits while handling the Boolean
identity/corner cases. For each segment it then:

1. allocates two four-integer vectors for input endpoint maxima;
2. forms eight 1,000-bit masked temporaries;
3. calls scalar `msb` eight times;
4. replaces absent classes with a large negative sentinel;
5. allocates a third four-integer vector for output extrema;
6. evaluates the endpoint-class maximum formulas;
7. fills parity-closed output sets one proxy bit at a time.

The eight searches at
[monitoring.hpp:4199](../../include/monitoring.hpp#L4199) dominate real random
inputs. With maximum input index 13, roughly 7,888 of their first 8,000 bit
tests per segment are guaranteed high zeros. The masked bitset operations
themselves touch 16 machine words each before those scalar loops begin.

The three `vector<int>(4)` allocations are needless, but replacing them by
`array<int, 4>` alone barely changes elapsed time because it does not remove
the reverse searches.

### Long-input timings and ablation

The following values are microseconds per segment. “Array” changes only the
three tiny vectors. “One scan” recovers all four endpoint extrema for a
language in a single reverse index pass, while retaining the public by-value
API, nested result, and scalar output fills.

| S | AND current | AND array | AND one scan | OR helper current | OR dual/move | OR dual+array | OR dual+one scan |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 64 | 4.135 | 4.123 | 2.536 | 4.243 | 4.217 | 4.186 | 2.542 |
| 128 | 4.195 | 4.167 | 2.532 | 4.236 | 4.214 | 4.196 | 2.555 |
| 256 | 4.124 | 4.047 | 2.613 | 4.146 | 4.079 | 4.069 | 2.626 |
| 512 | 4.173 | 4.110 | 2.601 | 4.202 | 4.109 | 4.099 | 2.600 |
| 1,024 | 4.147 | 4.125 | 2.646 | 4.306 | 4.228 | 4.199 | 2.553 |

Current conjunction is almost perfectly linear in `S`, at about 4.1
us/segment. Combining the eight searches into one gives 1.57--1.66x, but is
still scalar and still begins at index 999. A word-level implementation that
finds the highest nonzero 64-bit word and then uses leading-zero count should
be materially better, especially when an endpoint class is absent.

Output filling is not the dominant cost for the measured random languages:
AND results ended below index 25. It can matter when results approach capacity,
where parity masks and word-wise prefix fills should replace proxy-bit loops.

### Exact copy/allocation accounting

For normal lvalue calls on this ABI:

| Variant | Allocations per call | Requested bytes per call |
| --- | ---: | ---: |
| Current AND | `6S + 3` | `888S` |
| AND with arrays/one scan | `3S + 3` | `840S` |
| Standalone OR helper on lvalues | `8S + 5` | `1448S` |
| Expanded random-driver OR / helper with moved inner negations | `6S + 3` | `888S` |
| OR move plus arrays/one scan | `3S + 3` | `840S` |

At `S = 1,024` this was:

| Variant | Allocations | Requested bytes |
| --- | ---: | ---: |
| Current AND | 6,147 | 909,312 |
| AND array/one scan | 3,075 | 860,160 |
| Standalone OR helper on lvalues | 8,197 | 1,482,752 |
| Expanded random-driver OR / helper with moved negations | 6,147 | 909,312 |
| OR move plus array/one scan | 3,075 | 860,160 |

The `48S` bytes removed by arrays are transient integer blocks. The much
larger `840S` floor in the prototype is two copied inputs plus the nested
output. A const-reference kernel that handles index-zero cases without
clearing its inputs can avoid both input copies; with the current nested output
its allocation floor is `S + 1` and `280S` bytes.

The current identity unions at
[monitoring.hpp:4188](../../include/monitoring.hpp#L4188) use
`x = x | y`; `x |= y` avoids an explicit temporary at source level. This is
a small cleanup, not the main fix.

### Capacity bug

The fill loops at
[monitoring.hpp:4231](../../include/monitoring.hpp#L4231) do not check that a
computed index is less than `N`. For example, if both inputs contain bucket 1
at index 998, both `1...1` lengths are 999. The formula computes:

```text
len[3] = 999 + 999 - 1 = 1997
first written output index = len[3] - 1 = 1996
```

That calls unchecked `bitset::operator[](1996)` on `bitset<1000>`, which is
undefined behavior rather than a controlled capacity failure. The kernel
should validate every output extremum before filling, just as
[`bitsetConcat`](../../include/monitoring.hpp#L1693) checks its index.

### Recommended conjunction rewrite

In descending value order:

1. Take inputs by const reference and exclude singleton index zero in the
   extrema query rather than clearing copied inputs.
2. Replace eight masked scalar searches with word-level extraction of the four
   endpoint maxima. A cached four-extrema segment descriptor is even better if
   it can be maintained centrally.
3. Use `array<bitset<N>, 2>` for segments.
4. Fill parity-closed results by masked machine-word ranges.
5. Replace heap integer vectors with `array<int, 4>`, use `|=`, and add
   capacity and size checks.

The one-scan prototype matched current output at every tested size on the
repeated real languages. A production word-level rewrite still needs
exhaustive small-language and randomized differential tests.

## Disjunction

[`bitsetDisjunction`](../../include/monitoring.hpp#L4264) implements De Morgan:

```cpp
return bitsetNegation(
    bitsetConjunction(bitsetNegation(v1), bitsetNegation(v2)));
```

The random dispatcher expands that expression at its call site instead of
calling this helper. Starting from reusable lvalue APs, the expansion copies
each input once into its negation, then moves the two negated prvalues into
conjunction. Its OR shell therefore matches the measured “dual/move” variant:
`6S + 3` allocations and `888S` requested bytes, rather than the
standalone lvalue helper's `8S + 5` and `1448S`.

The mathematics is fine; the data movement is not. The public parameters are
already by-value copies. They become named lvalues, so each inner
`bitsetNegation(v1)` makes another deep language copy before swapping both
128-byte bitsets in every segment. The final negation performs a third full
swap pass over the conjunction output.

Changing the two inner calls to `bitsetNegation(std::move(v1))` and
`std::move(v2)` removes two deep copies. It reduced allocation traffic from
`1448S` to `888S`, but improved time by only about 1--2% in this workload:
the conjunction searches still dominate. This is a useful cleanup whose
importance grows after the search kernel is fixed.

The better design is one generalized binary kernel with a polarity flag or
bucket view. OR can then read each input with logical buckets reversed, run the
same endpoint algebra, and write reversed output buckets without materializing
three negated languages. That also permits const-reference inputs.

## Untimed always and eventually

### The central wasted computation

[`bitsetAlways`](../../include/monitoring.hpp#L2736) and
[`bitsetEventually`](../../include/monitoring.hpp#L2850) each allocate a complete
nested output eagerly. On every active-carry segment, they then allocate
`vector<int>(4)`, form four masked bitsets, and perform four scalar `msb`
searches.

Unlike conjunction, these functions do not need the maxima. Every subsequent
comparison asks only whether an endpoint class is empty, with one singleton
exclusion. The exact replacement is therefore a small set of presence
predicates.

Let `L0` and `L1` be the two input buckets, `even` and `odd` the parity
masks, and `nonSingletonEven` the even mask with bit zero cleared.

For always, while the suffix can still be true:

```text
out[0][0] = any(L0 & even) || any(L1 & odd)
out[0][1] = any(L0 & odd)  || any(L1 & nonSingletonEven)
out[1][0] = L1[0]
```

For eventually, while the suffix can still be false:

```text
out[0][0] = L0[0]
out[1][0] = any(L0 & odd) || any(L1 & even)
out[1][1] = any(L0 & nonSingletonEven) || any(L1 & odd)
```

These equations are simply the Boolean form of the existing length tests at
[monitoring.hpp:2758](../../include/monitoring.hpp#L2758) and
[monitoring.hpp:2872](../../include/monitoring.hpp#L2872). They retain the current
carry update exactly.

Do not reuse `allExceptFirstMask` without fixing initialization. The random
driver initializes only `evenMask` and `oddMask` at
[random.cpp:296](../../benchmarks/offline/approximate/random.cpp#L296);
`allExceptFirstMask` remains its zero-initialized global value. A rewrite
must explicitly initialize the required mask or construct a dedicated
`nonSingletonEven` mask.

### Carry-dependent behavior

For always, the incoming expensive-branch carry starts true and its next value
is `L1[0]`. Eventually is dual: its incoming carry starts true and its next
value is `L0[0]`. The first segment lacking the sustaining singleton is still
processed expensively; the branch then remains inactive for every earlier
segment. Thus `A` is the rightmost sustaining run plus its first failing
segment, or all `S` segments if no failure occurs.

The current work is therefore approximately:

```text
construct and zero the result for all S segments
+ 4A masked scalar reverse searches
+ A tiny heap allocations
```

This distinction matters. The tested random AP languages commonly produced
one of the two extremes: `A = S` or `A = 1`.

### Controlled results through 1,024 segments

Times below are microseconds per complete call. “Active-low” sustains the
carry and uses witnesses only at indices 0--14. “Collapsed” runs the expensive
branch once.

| S | Always active: current -> predicates | Always collapsed: current -> predicates | Eventually active: current -> predicates | Eventually collapsed: current -> predicates |
| ---: | ---: | ---: | ---: | ---: |
| 64 | 126.726 -> 1.505 (84.2x) | 3.313 -> 1.277 (2.59x) | 130.946 -> 1.701 (77.0x) | 3.326 -> 1.366 (2.43x) |
| 128 | 264.646 -> 3.347 (79.1x) | 4.681 -> 2.600 (1.80x) | 264.296 -> 3.620 (73.0x) | 4.611 -> 2.606 (1.77x) |
| 256 | 530.439 -> 6.353 (83.5x) | 7.460 -> 5.405 (1.38x) | 517.830 -> 7.083 (73.1x) | 7.558 -> 5.441 (1.39x) |
| 512 | 1,037.883 -> 13.141 (79.0x) | 13.276 -> 11.084 (1.20x) | 1,039.517 -> 14.643 (71.0x) | 13.433 -> 11.361 (1.18x) |
| 1,024 | 2,091.689 -> 25.044 (83.5x) | 22.737 -> 21.231 (1.07x) | 2,094.215 -> 27.408 (76.4x) | 23.127 -> 21.012 (1.10x) |

Active current cost remains about 1.98--2.07 us/segment; the predicates take
23.5--26.2 ns/segment for always and 26.6--28.6 ns/segment for eventually.
The stable ratio through 1,024 segments shows that this is not a short-input
artifact.

The inverse experiment confirms the cause. If all four endpoint classes have
a bit near indices 996--999, or the bitsets are dense, current `msb` stops
almost immediately. Current then costs about 32--37 ns/segment and the direct
predicates 24--29 ns/segment, only 1.25--1.55x faster.

### Actual long random input

A 1,024-sample, epsilon-8 random signal pair produced 968 canonical segments.
AP 0 supplies the active-always and collapsed-eventually rows; AP 1 supplies
the collapsed-always and active-eventually rows:

| Operator/carry state | Current | Direct predicates | Speedup |
| --- | ---: | ---: | ---: |
| Always, active | 1.9249 ms | 24.817 us | 77.6x |
| Eventually, active | 1.9269 ms | 26.384 us | 73.0x |
| Always, collapsed | 22.908 us | 20.385 us | 1.12x |
| Eventually, collapsed | 22.181 us | 20.360 us | 1.09x |

The active/collapsed split is more informative than a single average.

### Instruction and allocation evidence

Callgrind over five 512-segment active-low calls:

| Operator | Current instructions | Direct nested | Direct flat |
| --- | ---: | ---: | ---: |
| Always | 111.203 M | 1.580 M | 0.905 M |
| Eventually | 111.203 M | 1.718 M | 0.977 M |

That is 70.4x and 64.7x fewer instructions for the compatible nested-result
predicate implementations. For favorable high-index always input, the
instruction ratio falls to 1.41x, again isolating the reverse search as the
cause.

Allocation accounting is exact for nonempty input on this ABI:

```text
current:
  allocations = S + 1 + A
  requested bytes = 280S + 16A

direct predicates, nested result:
  allocations = S + 1
  requested bytes = 280S

direct predicates, vector<array<bitset<1000>, 2>> result:
  allocations = 1
  requested bytes = 256S
```

At `S = A = 512`, current makes 1,025 allocations and requests 151,552
bytes; the direct nested version makes 513 allocations and requests 143,360
bytes. Flat output makes one 131,072-byte allocation. Flat output was another
3.35--5.63x faster than direct nested output around this size, showing that
result construction is the floor after predicate work is removed.

There is an even deeper mismatch: each untimed unary result has at most three
meaningful bits per segment but materializes and zeroes 2,000 bits. A compact
tagged result, or a consumer that avoids materialization, would approach the
true lower bound. Flat generic storage preserves the encoding and was
semantically exact in the probe, but requires signature and caller changes.

### Validation of the direct predicates

The direct nested and flat implementations matched current always/eventually
exactly on:

- all 4,096 two-segment languages formed by every ordered pair of segment
  languages over both buckets at indices 0--2;
- 20,000 randomized languages of 1--128 segments, using 4 or 80 random bit
  draws per 2,000-bit segment across full-width indices 0--999;
- every controlled and real timing input above.

This covers singleton handling, every endpoint class, active and collapsed
carries, and full-width indices. The predicate rewrite is therefore the
best-supported optimization in this report.

## Timed always and eventually

### Current pipeline

The future bounded wrappers at
[monitoring.hpp:3123](../../include/monitoring.hpp#L3123) and unbounded wrappers
at [monitoring.hpp:3151](../../include/monitoring.hpp#L3151) take both the language
and segmentation by value, although
[`bitsetTimedUnary`](../../include/monitoring.hpp#L2906) reads them through const
references. A normal lvalue call therefore begins with a deep language copy
and a segmentation copy.

At the current random call site, the bounded-always wrapper instead receives a
negated language as an rvalue, so its language parameter moves. The deep copy
and bucket-swap pass have already occurred while constructing that negation;
the wrapper still copies segmentation. Unless stated otherwise, the isolated
public-wrapper tables below use an lvalue language and therefore include the
wrapper's language copy.

For every output segment, `bitsetTimedUnary` then:

1. scans every segmentation endpoint for every window offset to build
   `critical`;
2. creates singleton and open-interval placements;
3. constructs the timed window for each placement;
4. calls [`timedProfile`](../../include/monitoring.hpp#L1754), which scans every
   input segment, even those far outside the window;
5. concatenates every intersecting restricted segment into a complete
   alternating-word profile;
6. wraps that profile as a one-segment language and calls the full untimed
   always/eventually implementation merely to obtain two truth bits;
7. for ambiguous open placements, calls `changes` for every offset; each call
   scans all internal endpoints and all input segments again;
8. concatenates the placement result with generic `bitsetConcat`.

Always and eventually share this geometry. Their measured timed costs are
therefore almost identical.

### Bounds used by the random driver on longer inputs

The current random dispatcher uses future bounded always with `[0,1 s)` and
`[0,2 s)`; it does not directly call timed eventually. To isolate the unary
operator on its exact language, I constructed negated AP 1 before timing and
averaged five public lvalue-wrapper calls. Thus these values include the
wrapper copy, exclude construction of the negated AP, and are not complete
formula timings:

| Samples / canonical S | Always 1 s | Eventually 1 s | Always 2 s | Eventually 2 s |
| ---: | ---: | ---: | ---: | ---: |
| 64 / 60 | 3.016 ms | 2.963 ms | 5.416 ms | 5.591 ms |
| 128 / 127 | 8.417 ms | 8.710 ms | 17.143 ms | 16.843 ms |
| 256 / 249 | 17.796 ms | 16.242 ms | 33.364 ms | 31.610 ms |
| 512 / 477 | 25.784 ms | 24.965 ms | 54.916 ms | 54.892 ms |

The eventually columns show the shared operator geometry, not a call made by
the present dispatcher. The wider and unbounded experiments below are stress
tests used to expose the causes that these shorter bounds begin to exercise.

### Exact quadratic work for a fixed narrow window

With full-range evaluation over unit segments, an integer width
`1 <= w <= S`, bounded `[0,w)`, and distinct offsets `{0,w}`, every output
segment has exactly two placements: its left endpoint and its open interior.
Thus `C = 2S`.

The code performs exactly:

```text
critical candidate checks       = 2S(S + 1)
timedProfile segment checks      = 2S^2
included profile concatenations = (2w + 1)S - w^2
output-region concatenations     = 2S
```

The first two counts are quadratic and independent of window width. For
`w = 8`:

| S | Critical checks | Profile checks | Profile concats | Output concats | Total concat calls |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 64 | 8,320 | 8,192 | 1,024 | 128 | 1,152 |
| 128 | 33,024 | 32,768 | 2,112 | 256 | 2,368 |
| 256 | 131,584 | 131,072 | 4,288 | 512 | 4,800 |
| 512 | 525,312 | 524,288 | 8,640 | 1,024 | 9,664 |

If both truth values are possible on every open placement, `changes` is
called `2S` times. It adds:

```text
endpoint tests = 2S(S - 1)
segment tests  = 2S^2
```

The sum of these four loop-iteration counts is exactly `8S^2`: 2,097,152
iterations at `S = 512`, before profile and output concatenation work.

### Why a fixed window can look linear temporarily

The following synthetic stress timings for bounded `[0,8)` public wrappers
include the by-value language and segmentation copies. Sizes 64, 128, 256,
and 512 used
15, 8, 3, and 1 repetitions respectively, so each 512 row is a single
measurement:

| S | Constant-false always | Constant-false eventually | Ambiguous always | Ambiguous eventually |
| ---: | ---: | ---: | ---: | ---: |
| 64 | 3.050 | 3.097 | 11.890 | 11.664 |
| 128 | 6.513 | 6.449 | 24.746 | 24.610 |
| 256 | 13.101 | 12.753 | 48.358 | 48.157 |
| 512 | 27.766 | 27.671 | 99.485 | 99.695 |

These appear roughly linear through 512 despite the exact quadratic scan
counts. At these sizes, the `O(wS)` restricted-profile and concatenation work,
allocation overhead, and cache behavior likely dominate the cheap
outside-window intersection tests. The 3.6--3.9x ambiguous slowdown is
consistent with activation of the additional `changes` scans and `msb`
work.

The allocation traffic is already large. These `[0,8)` counts also include
the public wrapper's initial language and segmentation copies:

| S | Allocations/call | Requested bytes/call |
| ---: | ---: | ---: |
| 64 | 3,582 | 755,520 |
| 128 | 7,294 | 1,544,000 |
| 256 | 14,718 | 3,120,960 |
| 512 | 29,566 | 6,274,880 |

Counts are identical for always/eventually and definite/ambiguous inputs;
`changes` allocates no heap blocks. At `S = 512`, the 9,664
`bitsetConcat` results alone request 2,473,984 bytes, 39.4% of all requested
payload bytes. Allocator metadata is not included.

### The quadratic asymptote and width explosion are separate

A closed unbounded `[0,infinity)` window has a special untimed fast path at
[monitoring.hpp:3163](../../include/monitoring.hpp#L3163). Moving the lower bound
from zero to one selects the generic timed path and exposes clean quadratic
scaling on a constant-false unit language:

| S | Always `[0,inf)` | Eventually `[0,inf)` | Always `[1,inf)` | Eventually `[1,inf)` |
| ---: | ---: | ---: | ---: | ---: |
| 64 | 0.0105 ms | 0.152 ms | 10.60 ms | 10.45 ms |
| 128 | 0.0211 ms | 0.307 ms | 41.24 ms | 40.57 ms |
| 256 | 0.0598 ms | 0.633 ms | 161.67 ms | 161.35 ms |
| 512 | 0.186 ms | 1.320 ms | 647.46 ms | 648.71 ms |

The generic path grows by approximately 4x on every doubling. The fast-path
always/eventually difference is the untimed carry behavior described earlier.
The fast path itself still makes a second complete output and copies selected
segments into it; it can return the full result directly when the requested
range is the whole language.

Wider bounded windows expose the separate cost of actually concatenating
profiles. Each unit segment below contains only the constant-word language
`{0}` on even-numbered segments and `{1}` on odd-numbered segments:

| S | Width 1, A/E | Width 8, A/E | Width S/4, A/E |
| ---: | ---: | ---: | ---: |
| 64 | 1.009 / 1.012 ms | 3.117 / 3.128 ms | 5.280 / 5.172 ms |
| 128 | 2.119 / 2.108 ms | 6.494 / 6.565 ms | 19.769 / 19.689 ms |
| 256 | 4.904 / 4.933 ms | 14.063 / 13.776 ms | 75.532 / 75.293 ms |
| 512 | 12.594 / 12.692 ms | 30.996 / 30.392 ms | 297.506 / 297.731 ms |

For `w = S/4`, the number of included factor concatenations is
`(7/16)S^2 + S`, and the timing approaches 4x per doubling.

### Longer-window random-input stress

The real-input probe read `<samples>_0.txt` and `<samples>_100.txt`, used
duration `samples * 1000 ms`, epsilon `8000 ms`, and threshold 0, and passed
AP 0 alone to timed unary. All construction finished before the public-wrapper
stopwatch. The 64, 128, 256, and 512 rows used 5, 3, 1, and 1 repetitions:

| Samples | Canonical S | Window | Always | Eventually |
| ---: | ---: | ---: | ---: | ---: |
| 64 | 60 | `[0,8 s)` | 68.93 ms | 70.37 ms |
| 128 | 127 | `[0,8 s)` | 149.93 ms | 153.10 ms |
| 256 | 249 | `[0,8 s)` | 371.02 ms | 371.29 ms |
| 512 | 477 | `[0,8 s)` | 617.22 ms | 618.60 ms |
| 64 | 60 | `[0,16 s)` | 239.14 ms | 240.27 ms |
| 128 | 127 | `[0,32 s)` | 2,081.78 ms | 2,099.62 ms |
| 256 | 249 | `[0,64 s)` | 19,246.6 ms follow-up | 15 s run timeout |
| 512 | 477 | `[0,128 s)` | capacity exception | capacity exception |

The original 256-sample rows hit a 15-second process cap; a subsequent
single-call always run produced the value shown, while no exact eventually
value was retained. The 512-sample quarter-window run throws
`value expression exceeds the fixed bitset size` from generic profile
concatenation. This is a checked failure, unlike conjunction's unchecked
write, but it demonstrates that increasing `SIZE` would only postpone both a
time problem and a space problem.

Real AP factor languages and irregular segment boundaries make the timed path
much more expensive than the singleton synthetic case: they generate richer
intermediate profiles, more costly generic concatenations, and potentially
more placements. The timer still excludes all AP preprocessing.

### Specific timed bottlenecks

#### 1. Critical-point discovery

At [monitoring.hpp:2974](../../include/monitoring.hpp#L2974), every output segment
tests every endpoint against every offset, although a shifted endpoint can
belong to at most one half-open output segment. This is `Theta(kS^2)`
candidate work for `S` evaluated segments.

Bucket shifted endpoints into their containing output segment once, or merge
the sorted endpoint lists for each offset. A conservative first rewrite can
use `lower_bound` per segment while retaining `__int128` and exact
open/closed comparisons.

#### 2. Full segmentation scan for every profile

[`timedProfile`](../../include/monitoring.hpp#L1754) reconstructs a doubled
`TimedRange`, intersects it, and tests emptiness for every one of the `S`
segments on every call. A timed horizon is contiguous. Binary-search its first
and last intersecting segments and visit only that range. This removes the
`CS` outside-window checks while preserving the existing restriction
helpers.

That local fix does not solve wide or unbounded windows, because those really
touch many factors. The summary redesign below does.

#### 3. Constructing a full profile to discard almost all of it

At [monitoring.hpp:3048](../../include/monitoring.hpp#L3048), the code constructs
the entire alternating-word language of a window. It then allocates
`oneWindow`, calls a complete untimed operator, and keeps only two Booleans at
[monitoring.hpp:3070](../../include/monitoring.hpp#L3070).

For a profile `P`, the needed values are directly:

| Property | Exact test on a materialized profile |
| --- | --- |
| Can be all zero | `P[0][0]` |
| Can be all one | `P[1][0]` |
| Can contain zero | `P[0].any() || any(P[1] except index 0)` |
| Can contain one | `P[1].any() || any(P[0] except index 0)` |

Therefore:

```text
Always:     canFalse = canContainZero
            canTrue  = canBeAllOne

Eventually: canFalse = canBeAllZero
            canTrue  = canContainOne
```

More importantly, for the nonempty valid factor languages produced by the
monitor, these four properties compose over concatenation:

```text
allZero(L ++ R) = allZero(L) AND allZero(R)
allOne (L ++ R) = allOne (L) AND allOne (R)
hasZero(L ++ R) = hasZero(L) OR  hasZero(R)
hasOne (L ++ R) = hasOne (L) OR  hasOne (R)
```

So timed always/eventually need a four-Boolean monoid, not an alternating-word
profile. Compute summaries for the one or two restricted boundary factors and
use prefix counts for the full interior. Their summaries compose in constant
time after the boundary-factor summaries are available; the current
restriction helpers still inspect and copy `N`-bit storage. This removes
profile concatenations, most allocations, width-dependent language growth,
and the profile-concatenation `SIZE = 1000` overflow observed here. The
separate output-region bound can still exceed `N` and must retain its checked
failure.

#### 4. Repeated `changes` range scans

The lambda at [monitoring.hpp:2941](../../include/monitoring.hpp#L2941) scans all
internal endpoints, then all input segments, for every ambiguous placement and
offset. It also recomputes two scalar `msb` values for each intersecting
segment.

Precompute each segment's weight
`max(msb(bucket0), msb(bucket1))` once and its prefix sum. Count internal
endpoints using two bounds, locate the intersecting segment range, and answer
the weight sum as a range query. Saturate at `N` so the existing capacity
exception is preserved without integer overflow.

#### 5. Generic concatenation for highly structured truth regions

[`bitsetConcat`](../../include/monitoring.hpp#L1667) scans every index through
the left maximum. For every set left index and both right buckets, it recomputes
the same right `msb` at
[monitoring.hpp:1688](../../include/monitoring.hpp#L1688), then scans every right
index. Hoisting the two right maxima outside the left loop is an immediate
fix. Sparse-index lists or word-level shift/OR are larger improvements.

Timed unary output regions are much more structured than arbitrary languages:
constant zero, constant one, both singletons, or complete bounded prefixes.
They should have specialized composition, or be represented by endpoint
extrema, rather than sent through the generic Cartesian kernel.

After the four-Boolean window summary is implemented, generic concatenation
leaves the profile side but remains at
[monitoring.hpp:3113](../../include/monitoring.hpp#L3113) on timed unary's output
side until the specialized composition above is also implemented.

#### 6. Hot allocation and API details

Additional avoidable costs include:

- by-value public timed wrappers for read-only language and segmentation;
- a heap `vector<long long>` for at most two offsets;
- repeatedly growing `critical` and `placements` without useful reserves;
- a new profile, restricted factor, `oneWindow`, untimed result, region, and
  concatenation result for every placement;
- allocating a full `S`-segment output even when only `[s,e)` is evaluated;
- passing scalar `long long` and `bool` values by const reference.

Changing these details reduces traffic but does not cure the quadratic
algorithm.

### Timed complexity before and after

For `m` evaluated output segments, current geometry performs approximately:

```text
Theta(m k S)      critical candidate tests
+ Theta(C S)      full timedProfile scans
+ Theta(A k S)    ambiguous changes scans
+ profile and output concatenation work
```

Here `A` denotes ambiguous open placements in this formula, not the untimed
carry count. With `m = Theta(S)`, `C = O(kS)`, and fixed `k <= 2`, the
geometry is already quadratic. Wide/unbounded windows add
`Theta(S^2)` included-factor concatenations; a generic concatenation can
itself be quadratic in `N` in the dense worst case.

With shifted-event bucketing, four-Boolean range summaries, prefix-summed
`changes`, and specialized output composition, the expected structure is:

```text
O((S + C) * N/word)  whole-segment and boundary-factor summaries/weights
+ O(kS)              shifted critical events
+ O(C log S)         window/range lookup, or O(C) with monotone pointers
+ output materialization
```

That first bound assumes word-level scans; retaining the current scalar
restriction and `msb` helpers makes it `O((S + C)N)`. For fixed `k` and
`C = O(kS)`, the redesigned structure is near-linear in `S` and independent
of window width apart from output complexity.

This timed redesign was derived from the exact truth reductions but was not
implemented in this audit, so no speculative speedup number is assigned to it.
It requires careful differential testing of open/closed endpoints, punctual
windows, domain clipping, empty-window vacuity, restricted boundary segments,
and capacity behavior.

## Prioritized implementation plan

1. **Fix conjunction capacity checking.** This is correctness, not merely
   speed.
2. **Replace untimed always/eventually maxima with the validated predicates.**
   Initialize the non-singleton-even mask explicitly. This is small, local,
   and gives the clearest measured gain.
3. **Change fixed two-bucket segment storage to `array`.** It removes
   `S` heap allocations from every result and becomes important once scans
   are cheap.
4. **Rewrite timed unary around four Boolean window summaries.** Initially
   retain the existing exact restriction helpers and run old/new paths
   differentially before switching.
5. **Index timed ranges and `changes`.** Use boundary searches/prefix sums,
   then replace the per-output shifted-endpoint scan by a sweep or buckets.
6. **Make conjunction/disjunction const and word-oriented.** Use one generalized
   polarity-aware kernel; remove materialized De Morgan languages.
7. **Specialize timed output composition.** Only optimize generic
   `bitsetConcat` afterward if other operators still justify it.
8. **Apply low-risk cleanup:** arrays for tiny fixed data, moved negation
   arguments as an interim step, `|=`, reserves, and scalar-by-value
   parameters.

## Verification performed and remaining risk

The untimed direct-predicate variants received the strongest validation:
4,096 exhaustive two-segment low-bit languages and 20,000 randomized
full-width languages, for both always and eventually.

Binary one-scan/array/move variants were compared exactly with current output
at every controlled size on repeated real random languages. Broader exhaustive
testing is still warranted before changing the production binary kernel,
especially around absent endpoint classes, singleton identities, and the new
capacity check.

The repository's assertion-enabled timed regression executable passed during
this audit, but existing unary coverage is thin: it includes only one bounded
eventually equivalence and no comprehensive bounded-always or shifted
unbounded unary matrix. Before a timed summary rewrite, add exhaustive
differential tests against the old implementation for all four bound-closure
pairs, empty and punctual windows, clipped windows, and irregular
segmentations. Keep the old implementation available as a test oracle until
that matrix passes.

## Final assessment

The current speed on the original short random experiments hides three
different bottlenecks:

- binary operators do linear work with very large constants because low-index
  languages trigger almost full 1,000-bit scalar searches;
- untimed always/eventually calculate maxima they never use, producing an
  avoidable 70--84x active-path penalty;
- timed always/eventually repeatedly reconstruct full languages for questions
  that require only four Boolean summaries, combining quadratic scans with
  width-dependent profile growth.

Tiny-vector cleanup is worthwhile for allocation hygiene, but it is not where
most CPU time goes. The large improvements come from asking for exactly the
information each operator consumes: word-level endpoint extrema for AND/OR,
presence predicates for untimed always/eventually, and compositional Boolean
range summaries for timed always/eventually.
