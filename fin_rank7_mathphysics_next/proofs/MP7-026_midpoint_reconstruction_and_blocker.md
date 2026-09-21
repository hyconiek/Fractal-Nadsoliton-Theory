# MP7-026 — midpoint reconstruction, branch navigation, and exact proof blocker

Scientific state: **NUMERICAL_EVIDENCE / RECONSTRUCTION; MP7-026 REMAINS OPEN**.

The continuation checkpoint does not contain the source artifact named by MP7-025,
`R7P-031_simple_fold.json`, nor the original outward intervals for the four retained
spectral weights. Therefore a new proof-grade Lyapunov--Schmidt remainder bound cannot
honestly be produced from this checkpoint alone: higher-derivative enclosures would
otherwise depend on silently replacing certified parameter intervals by floating
midpoints.

A bounded reconstruction was nevertheless performed as a consistency check and to
reduce the remaining atom.

## Recovered midpoint model

Using the two separately certified nonzero stationary roots at exact `g=37/10` from
MP7-015 and the Fourier feature contract, a four-parameter nonlinear solve recovers

```
lambda3 = 1.96140686197644
lambda4 = 2.19956884933321
lambda5 = 2.29860627207909
lambda6 = 2.34218204114630
```

with maximum stationarity mismatch below `6e-16` at both roots. These values agree
with the midpoint values already recorded by the finite-N diagnostic artifact.

Solving the nine-equation stationary/null-vector/normalization system then recovers

```
s_fold = (1.36431142817824,
          1.43310270142855,
          1.40804570125504,
          1.00568207732017)

g_fold = 3.515644716839599

v = (0.507367539677250,
     0.528685655313811,
     0.553760694068750,
     0.395498105244207)
```

and the midpoint Hessian spectrum is approximately

```
(0, 0.17800274934, 0.21956286982, 0.23310794394).
```

The reconstructed fold coefficients are

```
a_mid = -0.212571633052236
b_mid =  0.118982597276860
```

which lie inside the already accepted MP7-025 intervals. Thus the reconstruction is
internally consistent with the checkpoint rather than a different local event.

## Numerical branch navigation

For `epsilon=g-g_fold` on a logarithmic grid `1e-7 <= epsilon <= 1e-2`, both local
stationary branches were solved directly. The dimensionless ratios were compared with
the MP7-025 leading coefficients:

```
xi_± / sqrt(epsilon)            -> ±1.890279088...
|Delta Phi| / epsilon^(3/2)     ->  0.535759617...
|lambda_soft| / sqrt(epsilon)   ->  0.224910315...
```

At `epsilon=1e-3`, the symmetric branch-separation ratio relative to the leading law is
within `~7e-5`, the local-energy ratio within `~2.5e-4`, and the two soft-curvature
magnitudes are within about 1.6% of the leading coefficient. This is a useful scale
choice for a future validated enclosure, but it is **not** the uniform remainder proof
requested by MP7-026.

## Smallest remaining proof atom

To promote MP7-026 to `PROVED_INTERVAL_ASSISTED`, import or re-supply the exact
R7P-031 fold box and the outward retained spectral intervals used by its checker. Then:

1. fix the normalized fold basis `(v,W)`;
2. use scaled coordinates `epsilon=t^2`, `s=s_fold+t z v+t^2 W w`;
3. interval-enclose the desingularized projected stationarity system through `t=0`;
4. prove a uniform inverse bound for its `(z,w)` Jacobian on `0<=t<=t_max`;
5. export explicit bounds for `z(t)`, `w(t)` and hence the three requested scaling errors.

The midpoint run shows that this is a local enclosure problem around a well-conditioned
three-dimensional transverse block, not evidence of a second hidden degeneracy.
