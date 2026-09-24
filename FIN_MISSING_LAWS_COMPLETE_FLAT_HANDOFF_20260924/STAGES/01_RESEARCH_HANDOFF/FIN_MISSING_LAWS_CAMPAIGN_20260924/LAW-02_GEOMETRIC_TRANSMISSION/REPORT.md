# LAW-02 — geometry-aware measure and transmission

**Disposition:** `PASS_ON_GEO002_CONTROL__GENERAL_RANDOM_GEOMETRY_OPEN`.

Take cell storage and edge transmission to be

`m_i = mu Vol(V_i)`,  `c_ij = kappa Area(F_ij)/ell_ij`,  `L=M^-1 K`.

On the exact rectangular-torus counterexample used by GEO-002, with `n_x=a m`, `n_y=m`, this gives

`L f = (kappa/mu) [delta_x^2 f / h_x^2 + delta_y^2 f / h_y^2]`.

Therefore both first-axis eigenvalues converge to `(kappa/mu) 4 pi^2`, independent of `a`. Numerically the weighted y/x ratio tends to 1 for `a=1,2,3,4`; the worst low-mode error has fitted slope about **-1.9983** in `m`, i.e. second order. By contrast the old unweighted density-scaled operator tends to ratio `a^2`.

This is more than “add ergodicity”: the *same nonergodic global aspect mixture* used in GEO-002 now has identical continuum operator in both components. Thus the earlier ergodicity/mixing obligation was specific to the unweighted RNG operator class.

The formula is the standard orthogonal finite-volume transmissibility `face measure / center distance`; it is externally known mathematics, not a uniquely FIN discovery. What remains novel for FIN would be sourcing the measure/dual-cell law from its relational state, and proving it on the actual state-generated geometry rather than a rectangular control.
