# R7P-073--080: cooperative structure with paid hypotheses

## Scope

The physical four-amplitude exponential family uses the finite abelian group
`G = Z4 x Z3` with uniform counting (Haar) measure and unscaled observables

- `c3 = cos(alpha)`, frequency `(1,0)`,
- `c4 = cos(beta)`, frequency `(0,-1)`,
- `c5 = cos(alpha-beta)`, frequency `(1,1)`,
- `c6 = cos(2 alpha)`, frequency `(2,0)`.

The physical Hamiltonian convention is `H=-sum_i J_i c_i`, with all `J_i>=0`.
This statement is about this shared-field model only.  The deliberately relaxed
split-`J5` negative controls are different Hamiltonians and are not covered.

## R7P-075: finite-group positive-definite proof of the hypotheses

For a character `chi` of a finite abelian group, let `v_g=chi(g)`.  Its
translation kernel is

`K_chi(g,h)=chi(g-h)=v_g conjugate(v_h)`,

hence is positive semidefinite.  Therefore

`K_Re(chi) = (K_chi + K_conjugate(chi))/2`

is positive semidefinite as well.  Each `c3,c4,c5,c6` is exactly such a real
part of a character.  Nonnegative linear combinations of positive-definite
functions remain positive-definite, so `-H=sum_i J_i c_i` is real
positive-definite whenever `J_i>=0`.

External ingredient: J. Ginibre, *General formulation of Griffiths'
inequalities*, Communications in Mathematical Physics **16** (1970), 310--328,
DOI `10.1007/BF01646537`.  Ginibre's general theorem together with Example 4
applies to the cone of real positive-definite functions on a compact abelian
group with Haar measure.  A modern finite-state restatement was also checked:
if `-H` belongs to this cone, both the first and second Griffiths inequalities
hold for cone observables.

All model-specific hypotheses have therefore been paid rather than inferred
from the word “ferromagnetic”.  Consequently, throughout the physical
nonnegative orthant,

`E[c_i] >= 0`,

and

`Cov(c_i,c_j) >= 0`

for all four observables and every pair.  Positive spectral normalizations
preserve all signs.  Thus for the fixed-point map

`T(s)=g E_s[C4]`

we have `T(s)>=0` and `DT(s)=g Cov_s(C4)>=0` entrywise.  Hence `T` is isotone
on the nonnegative amplitude orthant.

The earlier direct factorizations of `Cov(c3,c4)`, `Cov(c3,c6)`, and
`Cov(c4,c6)` remain independent exact checks.  The split-`J5` relaxed controls
remain useful demonstrations that the physical shared-field hypothesis matters.

## R7P-077: strict local Perron license

The universal theorem gives nonnegative entries, not strict positivity at every
boundary point.  On the validated localized and saddle root boxes, interval
evaluation is stronger: every covariance entry is strictly positive throughout
both boxes.  Therefore `DT` is a positive (primitive) matrix there and the
Perron eigenvector is simple and strictly positive.

The certified eigenvalue enclosures are recorded in
`certificates/R7P-077_local_jacobians.json`.  At the localized root all four
Jacobian eigenvalues are below one.  At the saddle exactly one is above one.
This is a statement about the declared fixed-point derivative, not an
unspecified physical dynamics.

## R7P-078: monotone fixed-point iteration

Let

`m=(sqrt(lambda3/6), sqrt(lambda4/6), sqrt(lambda5/6), sqrt(lambda6/12))`

be the coordinatewise maxima of the normalized `C4` features and set `B=g m`.
By R7P-075, `T` is isotone and nonnegative.  Also `T_i(s)<=g m_i=B_i` for all
`s>=0`.  Therefore `[0,B]` is invariant and `B` is a supersolution.
Uniform Fourier means vanish, so `T(0)=0`.

The iteration `s_(n+1)=T(s_n)` starting from `B` is componentwise decreasing
and bounded below by zero, hence converges coordinatewise.  Continuity of the
finite exponential family makes the limit a fixed point.  Any fixed point
`y in [0,B]` satisfies `y<=s_n` by induction, so the limit is the greatest fixed
point in the order interval.  Starting from zero remains at the minimal fixed
point zero.

This is an algorithmic convergence theorem for this fixed-point iteration only.
It is not a physical relaxation law.

## R7P-079: certified alignment diagnostic

At the exact rational gain `g=3.7183449`, both the localized and saddle roots
are interval-isolated.  The saddle Jacobian interval box supplies an operator
perturbation radius and a large certified leading spectral gap.  A Davis--Kahan
angle bound encloses the Perron direction, while the two Krawczyk boxes enclose
the localized-minus-saddle displacement.

The resulting cosine interval is stored in `certificates/R7P-079_alignment.json`:

`0.999554847191933 <= cos(angle) <= 0.9995548719371979`.

The midpoint angle is about `1.70963053` degrees.  This near alignment neither
proves an exact one-dimensional invariant path nor a minimum-action transition
trajectory.

## R7P-080 ledger conclusion

The original intake decision “universal cooperativity not accepted” was correct
at intake because the theorem hypotheses had not been checked.  They are now
checked.  The current ledger is:

- universal nonnegative means and covariance signs on the physical nonnegative
  orthant: **proved using Ginibre with explicit finite-group hypothesis check**;
- strict positivity / Perron simplicity at the two certified roots:
  **interval-certified locally**;
- monotone iteration from the explicit supersolution: **proved for the declared
  algorithm**;
- alignment: **certified scalar diagnostic, not a path theorem**;
- relaxed split-`J5` counterexamples: **negative controls outside the physical
  model**.

None of these conclusions implies the still-open global R7P-069 four-amplitude
curvature ceiling or resurrects the already-refuted everywhere full-seven-
coordinate index bound.
