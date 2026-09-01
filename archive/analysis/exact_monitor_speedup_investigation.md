# Exact-monitor speedup investigation

## Purpose

The RV-paper results reported a much larger speedup of the approximate
distributed monitor (ADM) over the exact SMT monitor (EDM) than the current
rerun. This note records the discrepancy, the experiments performed to explain
it, and the conclusions supported by the evidence as of 2026-08-22.

The largest discrepancy discussed here is for the random-trace formula
`always (p and q)` (`phi1`) at duration 32 and epsilon 1.

## Timing convention

An exact three-valued verdict requires running a positive and a negative exact
monitor instance. Its running time is therefore defined as

```text
max(positive-instance time, negative-instance time)
```

and not as the sum of the two times. Preprocessing is included in both recorded
instance times, following the exact benchmark driver. The historical results
also used the maximum rather than the sum. For example, historical sample 0 at
duration 32 and epsilon 1 records 14.5495 s for the positive instance and
13.7417 s for the negative instance; the relevant exact time is 14.5495 s, not
28.2912 s.

## Observed discrepancy

The corresponding heatmap-input cells are:

| Run | ADM | EDM | EDM / ADM |
|---|---:|---:|---:|
| Historical RV result | 0.00025149 s | 14.24021515 s | 56,623.39x |
| Current result | 0.00021172 s | 1.59534270 s | 7,535.15x |

The current cell was regenerated from the maintained benchmark outputs during
this investigation.

The change decomposes as follows:

```text
ADM improvement = 0.00025149 / 0.00021172 = 1.19x
EDM improvement = 14.24021515 / 1.59534270 = 8.93x
relative discrepancy = 8.93 / 1.19 = 7.51x
```

Thus, the approximate monitor became only about 19% faster, whereas the exact
monitor became almost nine times faster. The reduction from approximately
56,600x to 7,500x speedup is almost entirely due to the faster current exact
baseline.

The discrepancy is visible in the individual raw measurements and is not an
aggregation artifact. For sample 0, for example:

| Run | Positive exact | Negative exact | Verdict pair |
|---|---:|---:|---:|
| Historical | 14.5495 s | 13.7417 s | `(0, 1)` |
| Current | 1.4898 s | 1.1182 s | `(0, 1)` |

The current measurements are in the
[current reference](../../results/reference/offline/exact/random/phi1.txt)
raw result file.

## Implementations compared

The relevant sources are:

- Historical RV-paper implementation used for the recorded A/B test (not
  included here)
- [Current implementation](../../benchmarks/offline/exact/smt/random.py)

For the always-conjunction monitor and its dual negative check, the maintained
implementation reorganizes repeated code into helpers but retains the same
substantive discrete encoding: padding, signal constraints, inverse-clock
constraints, variability preservation, skew bounds, non-concurrent-edge
constraints, monotonic clocks, consistent-cut flow, and the final formula
query.

The current implementation additionally provides input validation, solver
timeouts, explicit handling of Z3 `unknown`, and clearer early-return logic. It
does not introduce solver reuse, incremental solving, parallel solving, or a
smaller encoding for this formula. It also omits the historical call to
`Solver.model()` after a satisfiable check.

## Controlled source A/B tests

To isolate source changes from hardware and solver-version changes, the
archived and current Python modules were run on this computer with the same
Python process environment, Z3 package, inputs, and timing procedure.

The comparison used:

- formula `phi1`, duration 32, epsilon 1;
- samples `0, 10, ..., 90`;
- a fresh child process for every implementation/sample measurement;
- one repetition per child;
- preprocessing included;
- `max(positive, negative)` as the exact time; and
- alternating archived/current execution order to reduce time-order bias.

### Z3 4.16.0

| Implementation | Mean exact time | Standard deviation |
|---|---:|---:|
| Archived RV | 1.60144 s | 0.21567 s |
| Current | 1.61377 s | 0.21095 s |

The archived/current ratio of means was 0.9924. The paired geometric-mean
ratio was 0.9918, with a 95% interval of 0.9788 to 1.0051. There were no
verdict mismatches among the ten inputs.

### Z3 4.13.0

The newest PyPI release available at the end of April 2024 was
`z3-solver==4.13.0.0`, released on 2024-03-07. The next release did not appear
until September. See the [PyPI release
history](https://pypi.org/project/z3-solver/#history).

It was installed without modifying the repository using a temporary `uv`
environment:

```bash
UV_CACHE_DIR=/tmp/z3-apr2024-uv-cache \
  uv run --with 'z3-solver==4.13.0.0' python ...
```

The environment resolved to Python 3.12.3 and Z3 4.13.0.

| Implementation | Mean exact time | Standard deviation |
|---|---:|---:|
| Archived RV | 1.69065 s | 0.24330 s |
| Current | 1.68013 s | 0.23810 s |

The archived/current ratio of means was 1.0063. The paired geometric-mean
ratio was 1.0056, with a 95% interval of 0.8962 to 1.1283. There were no
verdict mismatches among the ten inputs.

Under both tested solver versions, the archived and current source therefore
have indistinguishable performance. A refactor responsible for a ninefold
speedup would have produced a large difference in these same-machine A/B tests;
none was observed.

## Model extraction

The RV implementation calls `Solver.model()` after a satisfiable result,
whereas the current implementation only needs the satisfiability result. A
controlled measurement across two duration-32 samples and five repetitions
found that model extraction contributed approximately 0.000175 s to a total of
approximately 1.639 s, or about 0.0107%.

Omitting model extraction is a real cleanup but cannot explain the historical
speedup discrepancy.

## Repetition and warm-state effects

The effect of running the same monitor repeatedly in one process was also
measured. Across two duration-32 samples, the pooled three-repetition result was
only about 0.85% faster than the first run, with no monotonic warm-up trend.

Consequently, solver/Python warm state and the number of repetitions cannot
explain a factor close to nine.

## Z3-version effect

On the present machine, moving from Z3 4.13.0 to Z3 4.16.0 changed the mean by
only a few percent:

| Implementation | Z3 4.13.0 | Z3 4.16.0 | 4.16 improvement |
|---|---:|---:|---:|
| Archived RV | 1.69065 s | 1.60144 s | about 5.6% |
| Current | 1.68013 s | 1.61377 s | about 4.1% |

This test does not establish which Z3 version was installed on the RV-paper
computer. The paper states only that EDM invokes Z3 and contains no package
lock or exact version record. Version 4.13.0 is the correct latest-available
version for 2024-04-30, not proof of what was actually installed.

Nevertheless, Z3 4.13.0 runs both source implementations near the current
speed on the present computer. The 4.13-to-4.16 solver change does not reproduce
the old 14-second result.

## CPU and machine-environment assessment

The RV paper reports Ubuntu 24.04 on a laptop with an AMD Ryzen 7 4800HS and
16 GB RAM. The current machine reports an Intel Core Ultra 7 255U.

The processors are materially different:

- The [Ryzen 7 4800HS](https://www.amd.com/en/support/downloads/drivers.html/processors/ryzen/ryzen-4000-series/amd-ryzen-7-4800hs.html)
  is a Zen 2 mobile processor with boost frequency up to 4.2 GHz and 8 MB L3.
- The [Core Ultra 7 255U](https://www.intel.com/content/www/us/en/products/sku/241860/intel-core-ultra-7-processor-255u-12m-cache-up-to-5-20-ghz/specifications.html)
  is a newer mobile processor with performance-core boost frequency up to
  5.2 GHz and 12 MB cache.

Z3's configuration here is effectively single-threaded, so the comparison is
primarily about single-core execution rather than core count. As a rough scale
check, PassMark reports single-thread ratings around 2,545 for the
[4800HS](https://www.cpubenchmark.net/cpu.php?cpu=AMD+Ryzen+7+4800HS&id=3697)
and 3,649 for the
[255U](https://www.cpubenchmark.net/cpu.php?cpu=Intel+Core+Ultra+7+255U&id=6488),
a ratio of approximately 1.43. These third-party averages are not substitutes
for a controlled Z3 benchmark, but they show that an ordinary generational CPU
improvement is much smaller than 8.93x.

Z3 can benefit more than the approximate C++ monitor from improvements in
branch prediction, cache and memory latency, allocation performance, and
Python execution. The two workloads therefore need not scale identically.
However, explaining the observations by CPU alone would require the new CPU to
accelerate this Z3 workload about 7.5 times more than it accelerates ADM. That
is not credible as an ordinary model-to-model CPU difference.

An abnormal condition on the historical machine remains possible. Examples
include severe thermal or power throttling, virtualization or CPU quotas,
affinity to a slow core, heavy concurrent load, or a debug or otherwise poorly
optimized Z3 binary. Such conditions could affect a sustained, allocation- and
branch-heavy Z3 run much more than a roughly 0.2 ms C++ monitor invocation.

## Findings

The evidence supports the following conclusions:

1. **The timing aggregation is not responsible.** Historical exact timing used
   the maximum of the positive and negative instance times, not their sum.
2. **The current source refactor is not responsible.** Archived and current
   implementations run at essentially the same speed on the same computer and
   solver version.
3. **Model extraction is not responsible.** Its measured contribution is about
   0.01%.
4. **Warm-state and repetition effects are not responsible.** Their measured
   effect is below 1% in the tested configuration.
5. **The change from Z3 4.13 to 4.16 is not responsible.** It accounts for only
   about 4--6% on the present machine.
6. **A normal CPU-generation improvement is insufficient.** General
   single-thread performance suggests a much smaller difference, and ADM
   improved by only 1.19x.
7. **The discrepancy is specific to the historical exact-monitor execution
   environment.** The remaining candidates are the exact historical Z3
   version/build, Python/native-library stack, or an abnormal condition on the
   old machine.

The historical result need not be considered fabricated or arithmetically
wrong: it can accurately describe the old run while still being
machine/environment-specific and not reproducible on the current system.

## What would resolve the remaining uncertainty

The strongest follow-up would be to run the same archived source, inputs,
Python version, and pinned Z3 wheel on both the Ryzen 7 4800HS and the current
machine. For each machine, record:

- Python and `z3-solver` versions and package origin;
- the `libz3` binary hash and whether it is a release or debug build;
- CPU governor, power mode, temperature, and actual clock during the run;
- process affinity and virtualization/container limits;
- background CPU load; and
- separate expression-construction and `Solver.check()` times.

Hardware counters from `perf stat` (cycles, instructions, IPC, branches,
branch misses, cache misses, and page faults) would distinguish ordinary
single-core performance from throttling, memory behavior, or a substantially
different solver build. If the historical computer or its environment is no
longer available, testing Z3 4.12.6 and earlier releases on the current machine
can narrow the solver-version hypothesis, but it cannot reproduce or rule out
old-machine power, thermal, or software-build conditions.
