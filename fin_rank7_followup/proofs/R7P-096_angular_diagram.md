# R7P-096 locked angular fold/coexistence diagram

The constrained stationary system is formulated directly in orthonormal
Cartesian theta coordinates on the fixed-radius sphere, in the reflection-fixed
phase lock `(-pi/2,2pi/3,-pi/6)` with negative k6 sign.  This prevents mixing
the angular radius with the radial gain parameter.

## Full log-mgf
- fold: `r = 0.3463027188406942` (numerical bordered fold solve),
- equal-value event: `r = 0.3645555701282058`,
- pure-k6/localized equality is locally isolated by a six-dimensional interval
  Krawczyk box of radius `1e-8` at the fixed nominal strict spectral tuple,
- the intermediate locked angular saddle gives barrier
  `1.2902908116174672e-5`,
- the first pure-k6 instability is independently certified by R7P-095 near
  `0.41421132291`.

These reconstruct the imported full landmarks to the printed precision.

## Quartic truncation
The same constrained system for K4 gives
- coexistence `r = 0.3642826602833573`, consistent with the imported
  `0.36428264569` to about `1.5e-8`,
- fold `r = 0.346613544873972`.

The latter differs from the imported `0.34660000115` by about `1.35e-5`.
No historical generator for the imported number is present in the supplied
bundle, so the discrepancy is retained rather than tuned away.  It may reflect
an older solver/convention or a genuinely different quartic fixture and should
be resolved only from source provenance.

Scope: this is a locked angular diagram on a theta-radius sphere. It is not the
radial localization transition in g, a global rank-seven minimizer theorem, or a
physical-time hysteresis law.
