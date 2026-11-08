# Full nonlinear PDE parameter campaign

Date: 2026-11-08. Parent: 32bf806d42ea00e5a7225be64cb95317e60f09f5.

Thirty full nonlinear RP evolutions cover five w values, three amplitudes
and two coupled space/time resolutions. The fixed coefficients g4=g6=1,
physical preparation, T=epsilon^-2 and last-three-quarters inverse-action
slope estimator are unchanged from the pilot. No damping parameter is
fitted into the evolution. Canonical observables use the prior exact map.

The two resolutions are (dx,dt)=(0.2,0.04) and (0.1,0.02), with dt adjusted
to the endpoint. L=ceil(T+30). This campaign jointly refines space and time;
it does not independently separate those two errors at each new w. The
previous pilot's separate time/domain controls apply only to its w=1.5 case.

## Results on the finer mesh

| w | epsilon | Effective k_a | Relative difference from CD prediction |
| --- | --- | --- | --- |
| 1.25000000 | 0.20 | 0.416729890 | +1.7697% |
| 1.25000000 | 0.15 | 0.408846426 | -0.1555% |
| 1.25000000 | 0.10 | 0.408604635 | -0.2146% |
| 1.32333333 | 0.20 | 0.372137382 | +0.3948% |
| 1.32333333 | 0.15 | 0.368386724 | -0.6170% |
| 1.32333333 | 0.10 | 0.369799232 | -0.2360% |
| 1.33333333 | 0.20 | 0.366004920 | +0.1780% |
| 1.33333333 | 0.15 | 0.363052987 | -0.6299% |
| 1.33333333 | 0.10 | 0.364472320 | -0.2414% |
| 1.34333333 | 0.20 | 0.360026364 | -0.0094% |
| 1.34333333 | 0.15 | 0.357694754 | -0.6570% |
| 1.34333333 | 0.10 | 0.359145580 | -0.2540% |
| 1.75000000 | 0.20 | 0.190164026 | -1.5139% |
| 1.75000000 | 0.15 | 0.191458036 | -0.8437% |
| 1.75000000 | 0.10 | 0.192267535 | -0.4244% |

At epsilon=0.1, discrepancies are about 0.215% to 0.424%. Exactly at
w=4/3 the effective coefficient is 0.241% below the prediction. All runs
have positive continuum quadratic-energy gain. The maximum relative energy
drift on the finer mesh is 1.14e-7 (rounded upward). The largest paired
mesh change in k_eff/k_a is 0.000982 (rounded upward), or 0.0982 percentage
points. These are observed refinement changes, not certified error bars.

The amplitude trend is not monotonic in every case: at w=1.25 the
0.15-amplitude estimate is closer than the 0.1-amplitude estimate. All
results, including this feature and sign changes, are retained. No
asymptotic exponent is inferred from three amplitudes.

No abrupt threshold anomaly is resolved in this dominant finite-window
coefficient across w=4/3 +/-0.01. This does not rule out weaker cubic
threshold effects, longer-time changes or effects in the smaller shape
channel. This campaign principally probes k_a; it does not identify k_b
and k_x separately or validate a two-mode preparation theorem.

## Disposition and next step

The expanded campaign supports the leading transfer prediction in these
five cases. It is not independent-team replication, full-band validation,
a proof of AW or empirical FIN. AW remains conditional pending independent
review; routine self-audit chaining stays stopped.

Next implement a different numerical discretization and an explicit
interior outgoing-flux energy balance for a subset including the threshold
and both endpoints. Keep numerical uncertainty separate from finite-amplitude
and transient effects. Separately continue the strict-core nonlinear
matching problem already stated in FIN_BRIDGE_START.md; this campaign
supplies no derivation of the PDE from FIN.
