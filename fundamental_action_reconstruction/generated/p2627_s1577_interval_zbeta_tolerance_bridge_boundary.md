# P2627/S1577 interval Z_beta tolerance bridge boundary

Status: `P2627_INTERVAL_ZBETA_TOLERANCE_BOUNDARY_NO_STRICT_SOURCE_NO_FULL_BRIDGE_NO_LTOTAL_NO_QW2191_NO_TOE`

## Content-first anti-duplication audit

This packet greps research content rather than only ticket names/numbers before adding the tolerance result.

- `interval_bridge_content`: 20 hits
- `micro_coefficient_interval_content`: 6697 hits
- `nonlinear_damping_denominator_content`: 4492 hits
- `nonfit_noncircular_guard_content`: 1389 hits
- `closure_guard_content`: 13944 hits

## Theorem boundary

For the interval candidate D_Z(d)=1+(Z/100)d^(9/5) and strict denominator D_*(d)=1+d^(9/5), the relative error is |Z/100-1| d^(9/5)/(1+d^(9/5)).  On 0<d<=D this is maximized at D, so epsilon-admission requires |Z/100-1| <= epsilon*(1+D^(9/5))/D^(9/5).

For `epsilon=0.01` and `0 < d <= 10`, admitted `Z_beta` must lie in `[98.984151, 101.015849]`.
For a looser diagnostic `epsilon=0.15`, admitted `Z_beta` must lie in `[84.762266, 115.237734]`.

## QW-2064 coefficient rows

| source | Z_beta | sup relative error on d<=10 | pass 1% | pass 15% |
| --- | ---: | ---: | --- | --- |
| `target` | 100.000000 | 0.000000 | True | True |
| `micro_global_median` | 114.739580 | 0.145096 | False | True |
| `bin_q25` | 41.520217 | 0.575674 | False | False |
| `bin_q50` | 247.647250 | 1.453437 | False | False |
| `bin_q75` | 952.508925 | 8.392084 | False | False |

## Verdict

P2627 supplies the missing tolerance formula, but the available micro coefficient remains target-dependent and broad. At most it supports a future approximate/interval bridge lane; it does not export the exact positive_beta_renormalization_source required by P2625.

P2627 therefore does not repair P2625/P2620 and does not reopen role-bearing `L_total`, role transfer, `QW-2191`, or ToE closure.

## Recommended next honest step

Do not promote the interval lane to exact bridge closure.  Next, either narrow the micro Z_beta distribution with a target-independent operator/normalization law until the reported interval lies inside the chosen epsilon envelope, or explicitly downgrade the bridge program to an approximate effective-kernel theorem with declared epsilon and domain before any role-transfer audit.

Fingerprint: `a618a0825ec6a786f167b1636feed2f633db43305902f30b0d5788afb938af0f`
