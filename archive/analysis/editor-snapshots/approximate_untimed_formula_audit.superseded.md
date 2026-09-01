# Approximate Monitoring Efficiency Audit: Untimed Formulas

Date: 2026-08-22

## Scope and relation to the earlier audit

This report extends
[`approximate_random_operator_audit.md`](../approximate_random_operator_audit.md)
with the missing untimed random formula `phi3 = p until q`, a formula-level
comparison of `phi1` through `phi4`, and a word-oriented follow-up for
conjunction. Timed `phi5`/`phi6`, timed until, strict until, past operators,
preprocessing, exact monitoring, and the online architecture are out of
scope.

All formula timings below start after signal conversion, uncertainty
construction, canonical segmentation, value-expression construction, and AP
conversion. They measure only approximate formula evaluation. Microbenchmarks
used the same host and release-style `-O3 -DNDEBUG -march=native` setup as the
earlier audit. Comparisons in each row were compiled into the same binary.

No monitoring implementation was changed. Candidate kernels and allocation
probes were temporary analysis programs. The word-oriented prototypes use
libstdc++ set-bit iteration as a measurement proxy; production code needs an
explicit, supported word-level representation or accessor.

## Executive result

The two largest untimed opportunities are approximately tied in aggregate:

| Change | Long real-input result | Formula reach | Approximate saving |
| --- | ---: | --- | ---: |
| Const, cached, word-oriented conjunction kernel; direct polarity for OR | AND: 3.964 ms -> 156.3 us (25.4x) at 968 segments | `phi1`, `phi2`, `phi4` | About 3.8 ms per binary call, roughly 11.4 ms across one run of each affected formula |
| Const, cached, word-oriented non-strict until kernel | 11.746 ms -> 163.0 us (72.1x) at 968 segments | `phi3` | About 11.6 ms |
| Direct predicates for untimed always/eventually | Active always: 1.925 ms -> 24.8 us (77.6x) at 968 segments | outer/inner unary calls in `phi1`, `phi2`, `phi4` | About 1.9 ms per active-carry call |

For worst single-formula latency, until is the first target. For equal weight
over all four untimed formulas, the shared conjunction/disjunction kernel has
slightly broader impact. Both should use the same endpoint-summary primitive,
so the most useful first implementation step is shared.

Two correctness gates precede performance work:

1. add checked capacity failures to conjunction; and
2. add checked capacity failures to non-strict until.

Both functions currently reach unchecked `bitset::operator[]` with an index
outside the logical bitset capacity.

## Current untimed formula costs

The current random dispatcher evaluates:

| Formula | Current construction |
| --- | --- |
| `phi1` | `always(conjunction(p, q))` |
| `phi2` | `always(negation(conjunction(negation(p), negation(q))))` |
| `phi3` | `unboundedUntil(p, q, segmentation, 0, true)` |
| `phi4` | `always(negation(conjunction(p, always(negation(q)))))` |

Operator-only timings on the first real signal pair with epsilon 8 seconds
were:

| Samples / canonical segments | `phi1` | `phi2` | `phi3` | `phi4` |
| ---: | ---: | ---: | ---: | ---: |
| 64 / 60 | 0.254 ms | 0.364 ms | 0.709 ms | 0.362 ms |
| 128 / 127 | 0.516 ms | 0.762 ms | 1.504 ms | 0.766 ms |
| 256 / 249 | 1.010 ms | 1.022 ms | 2.956 ms | 1.990 ms |
| 512 / 477 | 1.934 ms | 2.872 ms | 5.756 ms | 2.880 ms |
| 1,024 / 968 | 3.932 ms | 5.829 ms | 11.560 ms | 5.876 ms |

At 968 segments, `phi3` is 42.5% of the total time for one evaluation of
each untimed formula. The approximately 3.9 ms `phi1` time is almost entirely
the current conjunction. `phi2` and `phi4` add an active untimed unary pass,
which explains their additional approximately 1.9 ms.

## Untimed non-strict until (`phi3`)

### Current call path and redundant copies

The dispatcher calls
[`bitsetUnboundedUntil`](../../../include/monitoring.hpp#L4027). Its `a = 0`, closed
lower-bound fast path:

1. copies both languages and the segmentation into by-value wrapper
   parameters;
2. calls [`bitsetUntilNonStrict`](../../../include/monitoring.hpp#L3638) with the two
   named lvalues, copying both complete languages again;
3. constructs a complete `full` result;
4. constructs a second complete result; and
5. copies the selected segments from `full` into that result.

The random `phi3` call selects the entire range, so the second result is
unnecessary. Segmentation is not used at all on this fast path.

The core takes its operands by value because it temporarily clears singleton
bits while handling identities. This mutation is avoidable: save the four
singleton flags as scalar Booleans and exclude index zero in the endpoint
query.

### Carry behavior and scan count

Non-strict until maintains two suffix carries, one for a false result and one
for a true result. Let `B0` and `B1` be the number of segments on which the
two branches run. The false branch starts active; after each segment, each
carry is simply whether the corresponding output bucket is nonempty.

On the real random inputs, both buckets remain possible almost everywhere:

| Samples / segments | `B0` | `B1` |
| ---: | ---: | ---: |
| 64 / 60 | 60 | 59 |
| 128 / 127 | 127 | 126 |
| 256 / 249 | 249 | 248 |
| 512 / 477 | 477 | 476 |
| 1,024 / 968 | 968 | 967 |

When the relevant singleton corner is present and the non-singleton algebra
also runs, each active branch performs:

- four scalar `msb` searches to reduce an always/eventually-like corner case;
- four scalar `msb` searches for the left operand; and
- four scalar `msb` searches for the right operand.

With both carries active, that is up to 24 reverse searches per segment. Each
search begins at index 999 although the real inputs end at indices 7--15.
The two branches recompute the same eight non-singleton endpoint extrema.

The function also allocates a `vector<bool>(4)` on every visited segment and
up to four `vector<int>(4)` objects per active branch. These are fixed-size
stack data disguised as heap data.

### Controlled scaling

The active-low input keeps both carries alive and uses only indices 0--14.
The word-oriented candidate replaces all repeated maxima with one cached
four-endpoint extraction per operand and segment. Times are per complete call:

| Segments | Current | Word-oriented candidate | Speedup |
| ---: | ---: | ---: | ---: |
| 64 | 766.4 us | 11.3 us | 67.7x |
| 128 | 1.535 ms | 22.6 us | 67.9x |
| 256 | 3.078 ms | 46.2 us | 66.6x |
| 512 | 6.156 ms | 91.1 us | 67.6x |
| 1,024 | 12.380 ms | 183.3 us | 67.5x |

The stable per-segment ratio shows that the gain is not a short-input
artifact. On deterministic singleton-only inputs, the current code reaches
the early `goto` without its main endpoint algebra. At 1,024 segments,
constant-false improved from 124.5 to 45.6 us (2.73x), and constant-true from
128.2 to 38.3 us (3.35x). This isolates the repeated reverse searches as the
source of the large active-low result.

### Real inputs

| Samples / canonical segments | Input max indices | Output max index | Current | Word-oriented candidate | Speedup |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 64 / 60 | 12 / 7 | 8 | 720.5 us | 7.3 us | 98.2x |
| 128 / 127 | 10 / 12 | 13 | 1.530 ms | 18.8 us | 81.5x |
| 256 / 249 | 13 / 12 | 13 | 3.011 ms | 42.0 us | 71.8x |
| 512 / 477 | 13 / 13 | 13 | 5.834 ms | 76.7 us | 76.1x |
| 1,024 / 968 | 14 / 15 | 15 | 11.746 ms | 163.0 us | 72.1x |

A portable scalar candidate that cached the same endpoint data but still
walked down from index 999 improved the 968-segment case only from 11.999 ms
to 2.613 ms (4.59x). The difference between 4.59x and 72.1x is direct evidence
that caching alone is insufficient: extraction must operate on machine words
or maintained metadata.

### Candidate structure

The validated candidate keeps the existing endpoint algebra and changes how
its inputs are obtained:

1. take both operands by const reference;
2. store the four singleton flags as Booleans instead of clearing input bits;
3. answer the always/eventually-like singleton corners with the direct
   presence predicates already validated in the earlier audit;
4. extract the four non-singleton endpoint maxima once per operand and
   segment;
5. cache those two summaries across both carry branches; and
6. replace all fixed-size `vector` objects with `array` or scalar fields.

The measurement prototype used libstdc++ `_Find_first`/`_Find_next` to
iterate set bits. Those names are non-standard and are not the recommended
public design. A production implementation should expose supported word
storage, scan the highest nonzero words with count-leading-zeros, or maintain
the four extrema as segment metadata.

### Allocation accounting

At 968 real segments, the direct current core and the const candidate made:

| Variant | Allocations | Requested bytes |
| --- | ---: | ---: |
| Current direct core on lvalues | 11,613 | 944,672 |
| Const cached candidate | 970 | 271,296 |
| Current public `phi3` wrapper | 14,522 | 1,765,800 |

The candidate's remaining allocations are almost entirely the nested output:
one inner allocation per segment. A `vector<array<bitset<N>,2>>` reduces that
floor to one allocation and improves locality.

The wrapper overhead is small while 24 scalar searches dominate, but it will
become a substantial fraction after the core rewrite. For the full-range
`phi3` fast path, return the core result directly and take the unused
segmentation by const reference (or bypass it entirely).

### Differential validation

Both the scalar-cached and word-oriented candidates matched the current
non-strict implementation on all in-capacity cases in:

- 16,384 exhaustive one-segment cases: every ordered operand language over
  both buckets at indices 0--2, under all four incoming carry pairs;
- 262,144 exhaustive two-segment cases: every operand sequence over both
  buckets at indices 0--1, under all four incoming carry pairs;
- 20,000 randomized languages of 1--64 segments, alternating between 4 and
  80 random bit draws per segment across the full 0--999 width; and
- the five real random inputs in the timing tables.

This covers empty classes, singleton identities, both carry branches, the
special sparse-copy cases, low and high indices, and dense inputs. It does not
validate out-of-capacity behavior because the current implementation has
undefined behavior there.

### Capacity defect

Non-strict until has an unchecked write distinct from the conjunction defect.
For one segment with incoming carries `(false, true)`, put:

```text
left bucket 1:  index 2
right bucket 1: index 999
```

Then `len1[3] = 3` and `len2[2] = 1000`. The weak-true formula computes
`len_weak[3] = len2[2] + 1 = 1001`, and the fill loop writes output index
1000 into `bitset<1000>` with unchecked `operator[]`.

On the tested libstdc++ layout, the invalid write lands in padding in the last
storage word: the resulting size-1,000 bitset reports `_Find_next(999) ==
1000` and `count() == 501`. This is still undefined behavior and violates the
logical encoding even when it does not cross an allocation boundary.
Production code must check every derived extremum before filling and throw a
controlled capacity exception.

## Conjunction word-oriented follow-up

The earlier audit measured a scalar one-pass conjunction prototype at only
1.57--1.66x because it still began at index 999. A const word-oriented
prototype, using the same endpoint algebra and no input mutation, produced:

| Samples / canonical segments | Current | Word-oriented candidate | Speedup |
| ---: | ---: | ---: | ---: |
| 64 / 60 | 244.7 us | 6.8 us | 36.0x |
| 128 / 127 | 524.6 us | 17.3 us | 30.3x |
| 256 / 249 | 1.035 ms | 39.7 us | 26.1x |
| 512 / 477 | 2.064 ms | 71.2 us | 29.0x |
| 1,024 / 968 | 3.964 ms | 156.3 us | 25.4x |

At 968 segments, allocation traffic fell from 5,811 allocations and 859,584
requested bytes to 969 allocations and 271,040 bytes. The remaining traffic
is the nested output.

The candidate matched current output on all five real inputs. Its algebra is
the already tested one-pass algebra from the earlier audit, but the new
const/word implementation has not received the same exhaustive randomized
matrix as the until candidate. Add that matrix before production use,
especially around absent endpoint classes, singleton identities, and capacity
failure.

Disjunction should not materialize De Morgan languages around this kernel.
Use a polarity-aware bucket view: read logically negated buckets, execute the
same endpoint algebra, and write logically negated output buckets. Once the
reverse searches are removed, the negation passes and deep copies become a
visible floor rather than a hidden 1--2% detail.

## Potential-impact order for untimed work

Correctness checks come first, but the performance work should be prioritized
as follows:

1. **Build one supported word-level endpoint-summary primitive, then apply it
   to a const, polarity-aware conjunction/disjunction kernel.** This affects
   three of four untimed random formulas. Direct AND is 25.4x faster on the
   longest real input and saves about 3.8 ms per call. Removing De Morgan
   materialization should improve `phi2`/implication further.
2. **Rewrite non-strict until around the same cached summaries.** This is the
   largest single-formula win: 72.1x and about 11.6 ms saved for `phi3` at 968
   segments. Replace corner maxima with predicates, make inputs const, and
   reuse each operand summary across both carries.
3. **Land the validated always/eventually predicate rewrite.** It saves about
   1.9 ms per active-carry pass on the long real input and affects `phi2` and
   `phi4` in the measured workload. Explicitly initialize or locally construct
   the non-singleton-even mask.
4. **Flatten segment storage to `array<bitset<N>,2>`.** After the computation
   fixes, nested output construction is the common floor: 969--970 allocations
   remain in the long conjunction/until candidates. A flat vector needs one.
5. **Remove full-range wrapper copies and duplicate materialization.** Return
   the full non-strict-until result directly for `phi3`; take unused
   segmentation and read-only operands by const reference. Apply equivalent
   ownership cleanup to negation and binary call paths.
6. **Consider compact/fused truth-shaped results.** Untimed unary outputs have
   at most a few meaningful bits per segment but zero 2,000 bits. A compact
   tagged result or direct consumer can remove the remaining materialization
   floor, but this is an architectural change and should follow the local
   exact rewrites.
7. **Apply small cleanup last:** fixed arrays, `|=`, reserves where relevant,
   and scalar parameters by value. These improve allocation hygiene but do
   not address the dominant scans.

Ranks 1 and 2 are a practical tie in total saving. If the goal is worst-case
latency for one formula, implement until first. If the goal is aggregate
throughput over `phi1`--`phi4`, implement the shared binary kernel first.

## Verification required before merging

For conjunction and until, add production differential tests before replacing
the current kernels:

- exhaustive small languages with every singleton combination and absent
  endpoint class;
- randomized sparse and dense full-width languages;
- every incoming carry pair for until;
- mismatched segment-count rejection;
- capacity-boundary cases that must throw before writing; and
- formula-level regressions for `phi1` through `phi4` on the fixed corpus.

Keep the current in-capacity implementation available as a temporary test
oracle. It must not be used as the oracle for capacity overflow, where its
behavior is undefined.
## Past-time port implementation results

The landed future-time implementation techniques have now also been applied
to untimed `Historically`, `Once`, and non-strict `Since`. This extends the
original audit's scope; it does not change the value-expression semantics.
Every result is still the complete pair of alternating-word bitsets, including
sparse non-singleton words, rather than only a pair of truth flags.

### Implemented architecture and invariants

- `bitsetAlwaysPast` and `bitsetEventuallyPast` are const-reference,
  left-to-right kernels. Direct bucket predicates replace repeated scalar
  `msb` searches and per-segment heap vectors. They emit the same singleton
  and length-two words as the old reductions.
- `reverseAlternatingSegment` reverses an alternating language without heap
  allocation. `reverseAlternatingSummary` preserves endpoint classes 0 and 3
  and swaps classes 1 and 2.
- The non-strict Until algebra is shared through one checked, allocation-free
  segment kernel. `bitsetSinceNonStrict` reverses each operand segment, calls
  that kernel, and reverses the complete output language back. Keeping the
  complete reversed bitsets is necessary for the Until kernel's sparse-copy
  branches; endpoint summaries alone would lose valid expressions.
- Closed `[0,infinity)` past wrappers evaluate the required prefix `[0,e)`
  before clearing `[0,s)`. This preserves history before `s`; evaluating only
  the requested slice with a fresh initial carry would be incorrect.
- Valid `bitset<1>` constant languages are supported, and every derived output
  length is checked before indexing. Malformed shapes and ranges now fail with
  deterministic `invalid_argument`, `out_of_range`, or `overflow_error`
  exceptions.

The carry transitions are unchanged. For `Historically`, the next false and
true possibilities are respectively `out0[0] || out1[1]` and `out1[0]`. For
`Once`, they are `out0[0]` and `out1[0] || out0[1]`. Non-strict `Since` keeps
the exact continuation rule

```text
nextFalse = sinceLastFalse || lhsLastFalse
nextTrue  = sinceLastTrue  && lhsLastTrue
```

The frozen hooks `bitsetAlwaysPastLegacy`, `bitsetEventuallyPastLegacy`, and
`bitsetSinceNonStrictLegacy` remain available as independent valid-input
differential oracles.

### Focused verification

The dedicated untimed-past target compares full languages, not only Boolean
verdicts. It covers:

- all 64 one-segment languages over indices 0--2, both unary operators, and
  all four incoming carries;
- all 16,384 ordered one-segment Since cases over the same language space and
  carries;
- all 1,024 languages over indices 0--4 for reversal involution and endpoint
  summary mapping;
- 240 randomized multi-segment operand pairs over every valid `[s,e)` and
  every carry, plus 500 randomized split/continuation cases; and
- history-sensitive partial wrappers, input immutability, `N = 1`, exact
  capacity, one-past-capacity rejection, and malformed calls.

The complete Release CTest run passed 7/7. The focused ASan/UBSan run passed
4/4 selected regression executables. Strict Since remains disabled under
`#if 0` and is absent from the maintained target and dispatcher lists.

### Measured performance

These are Release/O3 medians from seven samples after one warm-up. Each case
checked complete-language equality against its frozen legacy hook before it
was timed. `Ambiguous` means that every input factor contains both the false
singleton and the true singleton. `Rich` means both singletons plus one
length-two and one length-three alternating word; binary LHS and RHS fixtures
use opposite phases.

Representative full-range speedups at 256 and 512 segments were:

| Operator | Constant | Ambiguous | Rich |
| --- | ---: | ---: | ---: |
| Historically | 3.1x | about 105x | about 98x |
| Once | 1.3--1.9x | about 166x | about 165x |
| Non-strict Since | 2.43--2.87x across the measured patterns | 2.43--2.87x across the measured patterns | 2.43--2.87x across the measured patterns |

These are controlled synthetic operator measurements, not maintained
real-workload results.
