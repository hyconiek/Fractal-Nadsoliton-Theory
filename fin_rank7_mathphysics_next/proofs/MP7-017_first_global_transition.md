# MP7-017 — the certified equal-energy event is the first global transition

Scientific state: **PROVED_INTERVAL_ASSISTED_FIRST_GLOBAL_TRANSITION**.

## Statement

For each supplied strict spectral tuple covered by the accepted outward
intervals, let `(s_eq,g_eq)` be the unique local equal-energy stationary root
certified by R7P-026.  Then:

1. for every `0<g<g_eq`, the uniform state is the unique global minimizer;
2. at `g=g_eq`, the full global minimizer set consists of the uniform state and
   the twelve D12-related localized minima;
3. the localized branch has negative energy derivative through the event, so
   immediately above `g_eq` the uniform state is no longer globally minimizing.

Thus the R7P-026 event is not merely a local branch crossing: it is the **first
global coexistence/transition point** of the supplied finite rank-seven model.

The certified event enclosure is

`g_eq in [3.7183448971203875, 3.7183448991203876]`.

No physical provenance of g is asserted.

## 1. Global stationary exhaustion on the entire event parameter box

The global contractor of MP7-016 was rerun with g itself allowed to vary over
the full R7P-026 event interval.  In unscaled aligned coordinates

`J_i = g d_i m_i(J)`

with the accepted independent spectral enclosures, cooperativity again gives

`m_i(l) <= m_i(J) <= m_i(u)`

on every box.  Hence the contractor

`[l,u] -> [l,u] intersect [g_lo d_lo m(l), g_hi d_hi m(u)]`

is root-preserving simultaneously for every `(g,d)` in the event rectangle.

The outward-arithmetic tree has 253 records:

- 126 splits,
- 121 exclusions,
- one uniform terminal,
- one localized terminal,
- four terminal subboxes all contained in one common saddle branch tube.

Maximum depth is 76.  The complete tree is serialized in
`results/MP7-017_first_global_transition.json`.

On the uniform cube `[0,10^-8]^4` the interval infinity-norm row bounds for the
J-map derivative are at most

`(0.607765671831706, 0.681563048016401,
  0.712250990393494, 0.725753446913459)`.

Hence the uniform root is unique there.

## 2. Parametric stationary branch tubes

Separate four-dimensional Krawczyk tests, with g ranging over the entire event
box and the strict spectral tuple ranging over its accepted enclosure, give
uniform branch tubes for the localized minimum and saddle.

For the localized branch a radius-`10^-8` s-box strictly contains its Krawczyk
image.  The original 5D R7P-026 equal-energy box is contained in this branch
tube.

For the saddle, the corresponding radius-`10^-8` parametric Krawczyk box also
has strict inclusion.  Direct interval evaluation on the whole branch tube and
event g-box gives

`Phi_saddle >= 0.04655543192877835 > 0`.

Therefore, at the actual R7P-026 event of any supplied spectral tuple, every
aligned nonnegative stationary point is one of:

- uniform, energy exactly 0;
- localized, energy exactly 0 by the 5D stationary-plus-equal-energy
  certificate;
- saddle, energy strictly positive.

There is no fourth aligned stationary competitor.

## 3. Upgrade to full seven-dimensional global minimizers

MP7-011 proves that every full-X7 global minimizer is D12-equivalent to an
aligned nonnegative C4 representative.  MP7-014 excludes nonzero aligned
boundary stationary points throughout this gain window.  Hence every global
minimizer at `g_eq` must occur among the three exhausted aligned roots.

The saddle is excluded by its positive energy.  Thus the aligned global minima
are exactly the uniform root and the localized root.

MP7-034 proves that the localized aligned root has stabilizer of order two in
D12 and therefore exactly twelve distinct translated/reflected images.  The
uniform state is D12-invariant.  Consequently the full global minimizer set at
the event is

- one uniform state;
- twelve localized D12 images.

## 4. Why this is the *first* global transition

Use the primal potential

`V_g(p)=D(p||u0) - (g/2)||X7^T p||^2`.

At the certified event, globality just proved gives

`V_geq(p) >= 0`

for every simplex state p.  For any smaller gain,

`V_g(p)=V_geq(p) + (g_eq-g)||X7^T p||^2/2`.

Thus `V_g(p)>=0`.  If `p!=u0` and `X7^T p !=0`, the second term is strictly
positive.  If `X7^T p=0`, then `V_g(p)=D(p||u0)>0` unless `p=u0`.

Therefore the uniform state is the **unique** global minimizer for every
`0<g<g_eq`.

Finally R7P-028 certifies

`d DeltaV/dg in [-0.452413732176580,-0.452413730675052] < 0`

for the localized-minus-uniform branch difference at the event.  Hence for g
slightly above `g_eq` the localized branch has negative energy while the
uniform state remains at zero.  The event therefore genuinely changes the
global minimizing set.

## Nonconclusions

This theorem concerns the supplied dimensionless finite FIN model.  It does not
source a physical gain, temperature, clock, selector, laboratory transition,
SM/GR identification, or theory-of-everything closure.  It also does not claim
that the same localized orbit remains the unique nonuniform global phase for
all larger gains.
