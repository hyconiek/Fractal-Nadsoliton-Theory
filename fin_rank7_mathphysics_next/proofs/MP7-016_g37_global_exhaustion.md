# MP7-016 — complete aligned stationary exhaustion and global minimum at `g=37/10`

Scientific state: **PROVED_INTERVAL_ASSISTED_GLOBAL_EXHAUSTION**.

## Scope

This theorem is for the exact gain

`g = 37/10`

and the supplied strict spectral tuple enclosed by the accepted outward
intervals.  It uses the exact phase-alignment reduction of MP7-011 and the
boundary-support exclusion of MP7-014.

The complete machine-readable subdivision tree is
`results/MP7-016_g37_global_exhaustion.json` and is replayed by
`scripts/mp7_016_global_exhaustion.py`.

## 1. Unscaled-field stationary equation

Write

`J3=sqrt(lambda3/6) s3`, `J4=sqrt(lambda4/6) s4`,
`J5=sqrt(lambda5/6) s5`, `J6=sqrt(lambda6/12) s6`.

Let `c3,c4,c5,c6` be the corresponding unscaled cosine/alternating
observables and `m_i(J)=E_J[c_i]`.  The four stationary equations become

`J_i = g d_i m_i(J)`

with

`d=(lambda3/6,lambda4/6,lambda5/6,lambda6/12)`.

Since `0<=m_i<=1`, every aligned nonnegative root lies in the compact box

`0 <= J_i <= g d_i`.

Using the accepted spectral upper endpoints gives the explicit initial upper
corner

`(1.209534231552145, 1.356400790422153,
  1.417473867782112, 0.722172796020113)`.

## 2. Root-preserving isotone contractor

The accepted cooperativity theorem gives

`Cov(c_i,c_j)>=0`

throughout the nonnegative orthant.  Hence every mean `m_i(J)` is increasing in
every coordinate.  Therefore, on any box `[l,u]`,

`m_i(l) <= m_i(J) <= m_i(u)`.

For the independent accepted spectral enclosures this implies that any fixed
point in `[l,u]` must lie in

`[l,u] intersect [g d_lo m(l), g d_hi m(u)]`.

This intersection is a sound root-preserving contractor.  If it is empty, the
box contains no stationary root.  The checker evaluates the two corner means
with 60-decimal outward interval arithmetic, including the exact `sqrt(3)`
entries of the mode-five character table.  Every stored float endpoint is
inflated outward before interval evaluation.

Binary splitting is used only when the contractor does not yet decide a box.
A contracted-away shell is root-free by the displayed fixed-point implication,
so the tree does not need to cover that shell recursively.

The full run processes 253 records:

- 126 split records,
- 123 boxes excluded after contraction,
- 1 terminal box near the uniform root,
- 1 terminal box near the localized root,
- 2 terminal subboxes that both lie inside the same already-certified saddle
  Krawczyk box.

The maximum tree depth is 73.  Multiple terminal subboxes inside one Krawczyk
neighborhood do not represent multiple roots.

## 3. Local uniqueness of the three surviving neighborhoods

### Uniform root

On the explicit `J` cube `[0,10^-8]^4`, interval evaluation of the Jacobian of
`T(J)=g d m(J)` gives infinity-norm row-sum upper bounds

`(0.604767187226036, 0.678200475232203,
  0.708737014976567, 0.722172855512955)`.

Thus `||DT||_infinity < 1` throughout this cube.  Since `T(0)=0`, Banach
contraction gives the unique fixed point `J=0` there.

### Localized and saddle roots

For each surviving nonzero terminal box, the checker converts the entire
`J` box back to an outward `s` enclosure using every accepted feature scaling.
Each enclosure lies wholly inside the corresponding MP7-015 parametric
Krawczyk box.  Those source boxes already certify a unique stationary root for
the supplied spectral tuple.

Therefore the aligned nonnegative stationary problem has exactly three roots:

1. the uniform root;
2. one interior index-one saddle;
3. one interior localized local minimum.

There are no other aligned nonnegative stationary points at `g=37/10`.

## 4. Global energy ordering

Direct interval evaluation on the two nonzero Krawczyk boxes gives

`Phi_localized in [0.00823510200157903, 0.00823518137345460]`,

`Phi_saddle    in [0.0487900472912212, 0.0487901029016615]`.

The uniform root has exactly `Phi=0`.

MP7-011 proves that every global minimizer of the full seven-dimensional dual
is D12-equivalent to an aligned nonnegative C4 representative.  Such a global
minimizer is stationary and hence must be one of the three roots just
exhausted.  Both nonzero roots have strictly positive energy.

Therefore:

> At exact `g=37/10`, `theta=0` is the **unique global minimizer** of the full
> rank-seven dual.

By primal/dual minimizer correspondence, the uniform probability state is also
the unique global primal minimizer at this gain.

## 5. Immediate monotonic consequence

For fixed `p`,

`V_g(p)=D(p||u0) - (g/2) ||X7^T p||^2`

is nonincreasing in `g`.  Since `V_3.7(p)>=0` for every p, for every
`0<g<=3.7`

`V_g(p)=V_3.7(p)+(3.7-g)||X7^T p||^2/2 >=0`.

If `p!=u0` and the feature mean vanishes, the entropy term is strictly
positive; otherwise the added quadratic term is positive for `g<3.7`.
Consequently the uniform state remains the unique global minimizer for every
`0<g<=3.7`.

## Nonconclusions

This theorem alone does not yet prove that the certified local equal-energy
event near `3.71834489812` is globally first.  That requires the same global
exhaustion at the event itself, which is the next MP7-017 obligation.
