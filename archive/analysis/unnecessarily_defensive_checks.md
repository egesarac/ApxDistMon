# Review of Unnecessarily Defensive Checks in the Maintained C++ Entry Points

## Scope

There is no maintained literal `main.cpp`. This review covers the 17 executables
wired into the benchmark-root `CMakeLists.txt`:

- 7 incremental online drivers
- 7 naive online drivers
- 3 offline approximate drivers

The retired `archive/` tree is excluded. The disabled `since` drivers contain
the same online patterns, raising the source-copy count from 14 to 16.

## Summary

The maintained C++ entry points are mildly over-defensive, not severely so.
There are four families of provably dead or tautological checks, plus two
checks better expressed as debug assertions. More concerning is that online
input validation is repeated inside the timed region while still failing to
detect truncated input.

## Checks That Are Genuinely Unnecessary

### 1. `if (d < eps)` is unreachable

Representative instances:

- `benchmarks/online/incremental/main_eval_ac.cpp:197`
- `benchmarks/online/naive/main_eval_ac.cpp:196`

Every online driver restricts the parameters to:

- `d` in `{16, 32, 64, 128, 256, 512}` seconds
- `eps` in `{2, 4, 8}` seconds

Therefore `d >= 16 > 8 >= eps`; `d < eps` cannot happen. The check appears 14
times in built targets and 16 times including the disabled `since` drivers.

This check does not protect future arbitrary values either, because
`selectedMilliseconds()` rejects values outside those fixed sets before
reaching it.

**Verdict:** Remove it outright.

### 2. Several watermark checks are dead under the accepted configuration domain

The next filter requires:

```cpp
CHUNK_SZ > eps && CHUNK_SZ <= d
```

See `benchmarks/online/incremental/main_eval_ac.cpp:200`.

Consequently, for every configuration that actually runs:

1. On the first iteration, `dCurrent = CHUNK_SZ`.
2. Since `CHUNK_SZ > eps`, `qCurrent = dCurrent - eps > 0`.
3. Each later iteration strictly increases `dCurrent`.
4. Since `eps` is fixed, `qCurrent` also strictly increases.

This makes all of the following unreachable:

- The incremental backward-watermark exception at
  `benchmarks/online/incremental/main_eval_ac.cpp:266`.
- The incremental `qCurrent == qPrevious` branch at
  `benchmarks/online/incremental/main_eval_ac.cpp:280`.
- Its nested `qCurrent == 0` branch.
- The naive `qCurrent == 0` branch at
  `benchmarks/online/naive/main_eval_ac.cpp:233`.

There is a policy inconsistency here: the monitoring algorithm intentionally
supports chunks no larger than epsilon. Differential tests exercise
`chunk < eps` explicitly at `tests/test_online_monitoring.cpp:560`, while both
the C++ drivers and `scripts/run_online.py:207` filter such configurations out.

Two coherent choices exist:

- Preserve the paper's benchmark matrix and remove these dead runtime branches.
- Permit `CHUNK_SZ <= eps` in both the drivers and runner; then the zero/equal
  watermark branches become necessary.

The backward-watermark exception remains tautological either way:
`max(0, dCurrent - eps)` is monotonic whenever `dCurrent` is monotonic.

### 3. Refined-segmentation validation repeats a guaranteed helper postcondition

Representative incremental check:

```cpp
if (newSegmentation.size() < 2 ||
    newSegmentation.front() != qPrevious ||
    newSegmentation.back() != qCurrent)
```

See `benchmarks/online/incremental/main_eval_ac.cpp:316`. The equivalent naive
check is at `benchmarks/online/naive/main_eval_ac.cpp:252`.

The implementation of `refineSegmentation()` at `include/monitoring.hpp:1482`:

1. Inserts `left`.
2. Inserts only canonical/finalization endpoints strictly between `left` and
   `right`.
3. Inserts `right` whenever `left < right`.
4. Sorts and deduplicates.

Both callers invoke it only for a positive interval:

- Incremental evaluation already skipped the equal-watermark case.
- Naive evaluation already skipped the `qCurrent == 0` case.

It therefore necessarily returns at least `{left, right}`, with exactly those
endpoints at its front and back. The caller-side exception can fire only if the
helper violates its own construction.

**Verdict:** Remove it from Release benchmark code. If documenting the
invariant is desirable, use a debug `assert` or test `refineSegmentation()`
directly.

### 4. The naive final-endpoint check validates the same fact a second time

At `benchmarks/online/naive/main_eval_ac.cpp:265`:

```cpp
int last = root.size() - 1;
if (segmentation[last + 1] != qCurrent)
{
    throw logic_error(...);
}
```

For the current operators:

- `computeValueExpressions()` produces `segmentation.size() - 1` segments.
- Atomic propositions preserve that count.
- Every formula operator used by these drivers preserves segment count.
- The immediately preceding segmentation check already establishes
  `segmentation.back() == qCurrent`.

Therefore `last + 1 == segmentation.size() - 1`, making this another test of
`segmentation.back()`.

It is also not a robust defense against an operator regression:

- A shorter root is detected accidentally.
- A root longer than the segmentation can index past `segmentation.end()`
  before throwing.

**Verdict:** Remove it under current contracts. If structural validation is
wanted, the safe condition is:

```cpp
if (root.size() + 1 != segmentation.size())
```

That check belongs in a debug assertion or operator test.

### 5. The offline random formula fallback is unreachable

At `benchmarks/offline/approximate/random.cpp:218`:

```cpp
throw std::runtime_error("unsupported formula: " + formula);
```

`evaluateFormula()` has internal linkage and is called only at line 308, using
a `Formula` obtained from `selectFormulas()`. That selector can return only
`phi1` through `phi6`.

**Verdict:** Technically redundant. A better design would pass an enum instead
of a validated string. Retaining the fallback is harmless future-proofing and
is outside the timed region.

## Checks That Are Redundant but Useful as Debug Invariants

### Incremental cache alignment

At `benchmarks/online/incremental/main_eval_ac.cpp:358`, all formula caches are
compared with `retainedSegmentation.size() - 1`.

The invariant is inductively guaranteed:

- Initially there is one endpoint and zero cached segments.
- Each new segmentation contributes `k + 1` endpoints and each formula value
  contributes `k` segments.
- The code inserts `newSegmentation.begin() + 1`, contributing exactly `k`
  endpoints.
- Pruning erases the same `keep` prefix from every cache.

The `retainedSegmentation.back() != qCurrent` subcondition is especially
tautological because `newSegmentation.back()` was already guaranteed to equal
`qCurrent` and was just appended.

Still, this check catches future drift between formula operators and cache
management before `.back()` or pruning causes undefined behavior. Its O(1)
cost is tiny, but it occurs inside the benchmark timing.

**Verdict:** Keep it as a debug `assert` and compile it out of Release
measurements. At minimum, remove the redundant
`retainedSegmentation.back()` clause.

### Cross-repetition verdict consistency

At `benchmarks/online/incremental/main_eval_ac.cpp:420` and its naive
equivalent, later repetitions are compared with the first.

The monitor is deterministic, so this is theoretically unnecessary. In
practice it is valuable for catching hidden state, undefined behavior, and
accidental mutation. Importantly, the timer stops before the comparison, so it
does not contaminate reported runtime.

**Verdict:** Keep it.

## Misplaced Defense: Repeated Yet Insufficient Input Checking

The incremental driver checks `signalsReal[i].empty()` at
`benchmarks/online/incremental/main_eval_ac.cpp:229`. This is inside
`runMonitor()`, after the timer starts, so immutable input is revalidated during
every timed repetition.

That check is also insufficient:

- `getData()` at `include/monitoring.hpp:4631` silently returns empty data for a
  missing file.
- It returns short data for a truncated file.
- The incremental loop later indexes through `readEnd - 1` without checking the
  vector length.
- The naive drivers perform no precheck before `convertSignalsToBool()` at
  `include/monitoring.hpp:4668`, which also indexes up to the requested
  duration.

The offline random driver gets this right by validating exact lengths once,
outside the timed region, at
`benchmarks/offline/approximate/random.cpp:287`.

**Recommendation:** Validate both online input vectors once immediately after
loading, requiring exactly `d / 1000` samples, and remove the per-repetition
emptiness test.

## Checks That Should Remain

The following checks are not excessive:

- Missing option-value checks prevent `argv` overrun.
- `stoi` exception handling plus `parsed != value.size()` rejects both invalid
  integers and partially parsed strings such as `12junk`.
- Positive repetition checks prevent zero iterations and division by zero.
- Sample bounds protect the dataset naming convention.
- Supported duration/epsilon/chunk validation defines the benchmark matrix.
- Output-stream checks are needed because directory creation does not
  guarantee that file creation succeeds.
- The offline `epsilon > duration` filter is reachable for
  `duration=4, epsilon=8`.
- Offline exact signal-length validation prevents unchecked indexing.
- Agent-count validation is required before constructing `begin() + agents`.
- Mutual-separation coarse-variant rejection duplicates the lower-level
  evaluator, but the driver version fails before loading datasets or opening
  output files.
- `!retainedRoot.empty()` chooses between a temporal operator's initial state
  and carried state; it is semantic, not defensive.
- The top-level `catch (std::exception)` is appropriate at a CLI boundary.

## Recommended Cleanup Priority

1. Decide whether `chunk <= epsilon` is supported by benchmarks or only by the
   algorithms.
2. Remove the dead `d < eps` checks.
3. Eliminate the watermark branches made unreachable by the selected policy.
4. Move exact online input validation outside the timed lambda.
5. Remove caller-side segmentation postcondition checks.
6. Convert cache-alignment checks to debug assertions.
7. Replace string formula dispatch with a typed enum, making the unreachable
   fallback unnecessary by construction.
