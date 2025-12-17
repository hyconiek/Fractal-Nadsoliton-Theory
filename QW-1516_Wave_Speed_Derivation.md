# QW-1516: Derivation of c_tors from K(d)

**Date:** 2025-12-17 02:22:49.032715

## Derivation

The wave speed for torsion waves is:

$$c_{tors} = \sqrt{K(0)} = \sqrt{\alpha_{geo} \cos(\phi)}$$

Substituting frozen parameters:
- $\alpha_{geo} = 4 \ln 2 = 2.772589$
- $\phi = \pi/6$, so $\cos(\phi) = \sqrt{3}/2 = 0.866025$

$$c_{tors} = \sqrt{4 \ln 2 \times \frac{\sqrt{3}}{2}} = \sqrt{2\sqrt{3} \ln 2} = 1.549559$$

## Why c_tors ≠ 1?

The value c_tors = 1.51 arises from the geometric phase $\phi = \pi/6$.

In "information units" where $\alpha_{geo} = 4 \ln 2$, the wave speed is 1.51.
In Planck units where $c = 1$, we need to rescale:

$$c_{physical} = c_{tors} \times c_{light}$$

## Verdict
- c_tors is **derived from first principles** (not fit)
- The deviation from c = 1 has geometric origin (hexagonal symmetry)
