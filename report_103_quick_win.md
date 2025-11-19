# Report 103 — Quick-win RG probes
**Autor:** Krzysztof Żuchowski


Generated: 2025-11-14T17:29:18.621814Z

## Summary

This script computes simple RG proxies by scaling the distance kernel and measuring effective coupling, beta proxy, anomalous-dimension proxy, vacuum-energy proxy and related diagnostics. No fitting performed.

## Parameters

alpha_geo = 2.77

beta_tors = 0.01

omega = 0.7853981633974483


## Log

- TASK 0: N=16, alpha=2.77, beta=0.01, omega=0.7854

-   top eigenvalues (abs-sorted) = [20.751512, 11.028819, 11.028819, 0.682253, 0.682253, 0.036907, 0.036907, 0.012538]

- TASK 1-3: computed g(s) and numerical beta(g)

-   s_vals (sample) = [0.100, ..., 10.000]

-   g(s) [min,max] = [1.205976, 2.749458]

-   beta proxy (dg/dlns) [min,max] = [-2.258612e+00, 2.685457e+00]

-   near-zero beta indices = []

-   sign-change indices = [21, 26, 30, 32, 34, 35, 37]

- TASK 4: anomalous-dimension proxy gamma(s) estimated from top eigenvalue scaling

-   gamma(s) [min,max] = [-1.181083, 1.940470]

- TASK 5: vacuum-energy proxy E(s) computed (trace/N)

-   E(s) [min,max] = [2.770000e+00, 2.770000e+00]

- TASK 6: mass-hierarchy proxy (top/4th) across scales computed

-   ratio(s) [min,max] = [1.147460e+00, 1.083850e+04]

- TASK 7: operator-mixing proxy (overlap between top-3 eigenvectors at s=1 and s=2)

-   overlap matrix =
[[0.     0.     1.    ]
 [0.8409 0.5412 0.    ]
 [0.5412 0.8409 0.    ]]

- TASK 8: anomaly proxy (relative trace change s=1->2)

-   trA=4.432000e+01, trB=4.432000e+01, anomaly_proxy=0.000000e+00

- TASK 9: finite renormalization estimate (integrate beta between s=0.5 and s=2)

-   delta_g (approx) = -4.098767e-01

- TASK 10: phenomenological note on screening vs antiscreening from gamma(s)

-   screening_scales (sample up to 5) = [0.1, 0.1122, 0.1259, 0.1413, 0.1585]

-   antiscreening_scales (sample up to 5) = [1.2589, 1.4125, 1.5849, 1.7783, 3.1623]

