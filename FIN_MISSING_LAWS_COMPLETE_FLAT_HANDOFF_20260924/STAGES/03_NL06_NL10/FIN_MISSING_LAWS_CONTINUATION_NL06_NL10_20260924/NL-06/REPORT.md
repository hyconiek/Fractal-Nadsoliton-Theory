# NL-06 — actual FIN transverse potential on the current localized branch

Status: **PASS_UNIQUE_PHASE_QUOTIENT_MINIMUM_NO_TOPOLOGICAL_PHASE_SECTOR**

The previous NL-03 result was deliberately narrowed.  The full strict `T5`
phase torus has four invariant directions, but the currently accepted localized
continuum branch used here has only the active smooth harmonics `k=(3,4,5)`.
Hence its local phase quotient is two-dimensional.

Use

`alpha = phi4-phi3`,
`beta = phi3+phi5-2 phi4`,
`gamma = 4 phi3-3 phi4`.

Then

`phi3=3 alpha+gamma`,
`phi4=4 alpha+gamma`,
`phi5=5 alpha+beta+gamma`.

The all-orders positive-resonance theorem says that for positive amplitudes the
global equality manifold is exactly `phi_k=k alpha`.  Therefore, after quotient
by alpha, `(beta,gamma)=(0,0)` is the **only global minimum orbit** of the
phase potential.  This is an exact structural conclusion, not a grid claim.

At the accepted `g=7` localized amplitudes, the quotient Hessian is

```
 0.812795448174  -0.300518531529
-0.300518531529   0.179540089893
```

with eigenvalues `0.059630857661` and `0.932704680406`.
Both transverse directions are gapped.

A 300-start whole-torus numerical search found only the locked quotient
minimum.  Reflection gives `U(-beta,-gamma)=U(beta,gamma)`, so every odd total
Taylor order vanishes at the origin.  The quartic tensor is nontrivial (see
`results.json`), but it does not create a second degenerate well; higher orders
restore the globally unique quotient minimum.

**Consequence:** the current FIN phase-deformation sector does not source a
sine-Gordon-like pair of inequivalent vacua.  Searching for a topological kink
in `(beta,gamma)` would therefore be target fitting.  This conclusion is only
about the phase quotient; it does **not** exclude the distinct uniform/localized
radial phases studied later in NL-10.
