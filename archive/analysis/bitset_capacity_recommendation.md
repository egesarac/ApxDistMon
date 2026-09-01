# Bitset-capacity recommendation

## Decision

Computing and recording safe capacity bounds is worthwhile for correctness.
Building a complicated per-input capacity system is not currently worth the
maintenance cost.

The practical design should be:

1. Make the compile-time capacity configurable through
   `STTT_BITSET_SIZE`.
2. Keep `1000` as the default for general programs and tests.
3. Assign a smaller capacity to experiment executables that are already
   separate.
4. Keep the combined offline-random executable at `340` instead of splitting
   it into six executables solely for capacity tuning.
5. Add explicit overflow checks. A capacity error must produce a clear failure,
   not an unchecked `bitset::operator[]` write.
6. Reconsider per-formula offline executables only after the direct-predicate,
   word-level-extrema, and copy-removal optimizations have been benchmarked.

`std::bitset<N>` has a compile-time capacity. Consequently, a command-line
`--size` option is not useful unless the program contains several separately
instantiated monitor types and dispatches between them. That complexity is not
justified here.

## Recommended capacities

`SIZE` is a capacity, not a maximum index. If index `k` can be written, the
required capacity is `k + 1`.

### Random formulas

For the maintained random workloads, including durations longer than the
offline grid and clock skew up to eight seconds, the formula-specific bounds
are:

| Formula family | Minimum safe capacity |
| --- | ---: |
| `phi1`: `G(p and q)` | 33 |
| `phi2`: `G(p or q)` | 33 |
| `phi3`: `p U q` | 18 |
| `phi4`: `G(p -> F q)` | 18 |
| `phi5`: `G(p -> F_[0,1) q)` | 279 |
| `phi6`: `G(p -> F_[0,2) q)` | 340 |

The corresponding online targets use the same operator bounds:

| Online target | Minimum safe capacity |
| --- | ---: |
| `ac`, `ad`, `ec`, `ed` | 33 |
| `response`, `since` | 18 |
| `response01` | 279 |
| `response02` | 340 |

Both online modes use the same capacity for a given formula. Chunk size changes
which part of the trace is processed at once, but it does not increase the
largest local word that the operator can construct.

The offline random driver currently evaluates all six formulas from one
executable. Giving that executable capacity `340` is therefore the simplest
safe choice. Separate capacities of 33, 18, 279, and 340 would require either
several formula-locked executables or a larger type-dispatch refactor.

### Case studies

The case-study maxima arise before the final untimed `Always`. Numeric signal
languages are combined and then converted to a Boolean predicate. The returned
invariant is small, but the predicate must fit in the bitset first.

For the maintained experiment grids:

| Executable | Minimum safe target capacity |
| --- | ---: |
| Water tanks | 65 |
| Mutual separation | 61 |

More detailed water-tank bounds are:

| Variant class | Maximum safe capacity |
| --- | ---: |
| `adm`, `adm-f` | 65 |
| `adm-fr` | 49 |
| `adm-c`, `adm-cr` | 2 |

The single water-tank executable supports all of these variants, so it should
use `65`.

For mutual separation, `adm` and `adm-f` require at most `61`. The relative
variant needs only `31` for two agents, but it again reaches `61` with three or
four agents because there is then a pair containing two non-reference agents.
The shared executable should therefore use `61`.

These case-study bounds assume the maintained sampling structure and parameter
grid. A replacement dataset accepted through `--data-dir` could have more
changes. The executable must validate the input-derived requirement against
its compiled capacity and reject an unsupported dataset cleanly.

## Memory effect

On the measured GCC/libstdc++ platform, capacities are rounded to 64-bit
storage words:

| Capacity | Bytes per `bitset` | Reduction from `bitset<1000>` |
| ---: | ---: | ---: |
| 1000 | 128 | baseline |
| 340 | 48 | 2.67x smaller |
| 65 | 16 | 8x smaller |
| 61 | 8 | 16x smaller |
| 33 | 8 | 16x smaller |
| 18 | 8 | 16x smaller |

Each Boolean segment contains two bitsets, and the current operators copy and
allocate many complete languages. The payload reduction is therefore real,
although the number of heap allocations does not change.

## Expected speed effect

The following figures are expectations, not capacity-only benchmark results.
They should not be reported as measured speedups.

With the current scalar `msb` implementation, a low-index word starts its
search at index 999. On the untimed random workload, useful bits end below
approximately index 33. Reducing the capacity to 33 or 18 removes roughly
950--980 guaranteed-empty bit tests from each search. A 10--40x improvement in
the affected current operator kernels is plausible. The complete experiment
would improve less—roughly 2--8x is a reasonable range to test—because signal
loading, uncertainty construction, segmentation, value-expression creation,
and allocation remain.

Reducing the combined random executable from 1000 to 340 is much less
dramatic. It reduces bitset payload by 2.67x and shortens the current reverse
search by about threefold. A low-single-digit improvement is the realistic
expectation.

Capacity reduction does not compose multiplicatively with the main planned
operator improvements:

- Direct `Always`/`Eventually` predicates already remove the unnecessary
  maxima searches. The audit measured 70--84x on active-carry paths at
  capacity 1000.
- Word-level conjunction/disjunction extrema replace hundreds of scalar bit
  tests with a few machine-word operations.
- Removing input copies and materialized negations reduces the amount of
  bitset data moved.

After those changes, a formula-specific capacity is expected to add only a
small-to-moderate gain—probably about 1.2--3x in the affected operator, and
less in the complete experiment. This is why splitting the offline random
driver is premature.

For timed formulas, changing 1000 to 340 may improve the low-level bitset work
by roughly 1.5--3x, but it does not change the timed algorithm's quadratic
segmentation scans or its allocation-heavy profile construction. It is not a
substitute for redesigning timed unary evaluation.

For the case studies, the likely end-to-end gain is small. Their expensive
work is numeric string-language construction; the final Boolean predicate and
untimed invariant are a comparatively small tail. Smaller bitsets are still
useful for correctness, memory, and consistency, but should not be presented
as a major case-study optimization.

## Suggested implementation boundary

The worthwhile, low-complexity change is:

```cpp
#ifndef STTT_BITSET_SIZE
#define STTT_BITSET_SIZE 1000
#endif

inline constexpr std::size_t SIZE = STTT_BITSET_SIZE;
```

CMake can then apply `STTT_BITSET_SIZE=<capacity>` privately to each existing
experiment target. Tests and unclassified executables retain the default.

The implementation should also:

- check conjunction/disjunction output extrema before writing bits;
- check the requested experiment configuration or input-derived bound against
  the compiled capacity;
- include the compiled capacity in `--list` or diagnostic output;
- keep capacity out of timing-result rows unless the result-file schema is
  deliberately versioned;
- add boundary tests whose largest valid index is exactly `SIZE - 1`, followed
  by an input that must be rejected.

## When formula-specific offline binaries become worthwhile

Benchmark after the direct predicates, faster extrema, and copy-removal work.
Only split the combined offline random driver if using capacity 33 or 18
instead of 340 still produces a meaningful end-to-end improvement. A sensible
threshold is approximately 15--20%; below that, the additional targets,
formula locking, runner mapping, build time, and maintenance are unlikely to
pay for themselves.

The detailed operator evidence and measured microbenchmarks are in
[approximate_random_operator_audit.md](approximate_random_operator_audit.md).
