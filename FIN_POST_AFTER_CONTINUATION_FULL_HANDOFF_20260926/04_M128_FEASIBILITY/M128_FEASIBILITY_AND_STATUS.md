# GEOMETRY-M128 / ASYMPTOTIC DEPTH — FEASIBILITY AND CURRENT STOP

This continuation did real implementation/feasibility work toward the P0
`GEOMETRY-M128 / ASYMPTOTIC DEPTH` target, but it did **not** produce a
completed scientific M=128 result.  Do not promote M128 as closed.

## Durable implementation state

Portable sources retained in this handoff:
- `m128_q1_7d.cpp`
- `m128_q1_7d_float.cpp`
- `bench17_7_fftw.cpp`
- `bench17_7_fftwf.cpp`
- `make_q2_lane.py`

The 7-D FFT grid used by the producer has `17^7 = 410,338,673` grid points.
The in-place real FFT allocation is `24,137,569 * 18 = 434,476,242` scalars.

Observed benchmark:
- float FFT allocation: about **1.738 GB**, FFT about **3.75 s**;
- double FFT allocation: about **3.476 GB**, FFT about **5.06 s**.

The producer therefore became computationally feasible at the individual-FFT
level, but the full exact M128 campaign requires several such work arrays plus
parent/derivative state and was not completed in this continuation.

## Completed M64/depth replays retained

Representative float runs:

| alpha | C_H | root fraction | root deficit | internal C_H |
|---:|---:|---:|---:|---:|
| 0.850 | 373.4121664 | 0.61130788 | 0.38869212 | 145.1423672 |
| 0.875 | 992.1824249 | 0.86977168 | 0.13022832 | 129.2102500 |
| 0.950 | 112.3951877 | 0.86284324 | 0.13715676 | 15.4157594 |

These values do **not** justify a monotone root-deficit exponent from M<=64.
The original falsification criterion remains active.

Compact q2 metadata for alpha 0.85, 0.875 and 0.95 is included under
`03_COMPACT_DATA/`.

## Raw intermediates intentionally excluded

The working directory contains multi-gigabyte raw convolution arrays
(including ~1.7 GB M64 work files and many 19–37 MB derivative arrays).
They are not included in the handoff ZIP because they are regenerable
intermediates and would make the handoff multi-gigabyte.

Their presence is recorded in `EXCLUDED_RAW_INTERMEDIATES.csv`.

## Status

**OPEN / PARTIAL COMPUTATIONAL PROGRESS.**

Next valid action:
finish a resource-aware M128 scientific replay or formally declare a resource
stop.  No exponent should be promoted before that.
