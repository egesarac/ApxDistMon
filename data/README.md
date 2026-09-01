# Benchmark input data

These files are fixed benchmark fixtures. They allow exact replay of the
paper's monitoring inputs, but the original data-generation process is not
fully reproducible.

## Random-generator traces

`signals/` contains 1,800 files: 200 files for every duration
`d={4,8,16,32,64,128,256,512,1024}`. A file named `d_i.txt` has exactly `d`
whitespace-separated rows:

```text
time_seconds integer_value
```

Times are `0,...,d-1`; values are integers in `[-100,100]`. Benchmark sample
`c`, for `c=0,...,99`, uses `d_c.txt` and `d_(c+100).txt` as its two signal
components. Values are mapped to Boolean propositions with `value > 0`; zero
is false. The runners apply clock skew `epsilon` to these fixed traces, so
`epsilon` is not an input-generation parameter.

The generator used C++ `uniform_int_distribution(-100,100)` and seeded a new
`mt19937` from `random_device` for every file. It did not accept or record a
seed or engine state. The supplied files can therefore be replayed exactly,
but they cannot be regenerated bit-for-bit from documented parameters.

## Case-study traces

`case_studies/` contains the one fixed 1 s trace set used in the paper:

- `uav/s1_uav_i`, for `i=0,...,3`, has 21 rows
  `time_seconds x y z` at 20 Hz. The row at `t=1` is an endpoint marker, so the
  semantic trace has 20 intervals over `[0,1)`.
- `water_tanks/s5_tank_i`, for `i=0,...,3`, has 20 rows
  `time_seconds outflow_pressure` at 20 Hz over `t=0,...,0.95`. The final
  value is held through `t=1`.

A configuration with `n=2,3,4` agents uses files `0,...,n-1`. Coordinate and
pressure units are not documented.

According to the paper, the water-tank traces were exported from its Simulink
model and the UAV traces from the Fly-by-Logic path planner. The package does
not contain the source models, scenarios, generation/export scripts, software
versions, run seeds, or settings. These exported traces are replayable, but
their generation is not independently reproducible from this repository.
