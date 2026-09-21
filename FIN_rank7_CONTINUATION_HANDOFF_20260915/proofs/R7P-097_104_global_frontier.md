# R7P-097--104 — full-seven-coordinate global frontier

## Exact variational target

For nonuniform simplex probabilities with
`Q=(p-u)^T A7 (p-u)>0`, define

`R(p)=2 D(p||u)/Q`.

Then `V_g(p)=D(p||u)-g Q/2`, so a negative-energy competitor exists exactly
when `g>R(p)`.  The first energetic threshold is therefore the infimum of
`R(p)` over `Q>0`; `Q=0` states cannot lower the energy because their quadratic
term vanishes.  The tangent limit at uniform is a separate local quantity and
is not the finite-amplitude crossing.

## Certified initial bracket

The current repository guardrail ST448 proves that for the full strict
Laplacian the uniform probability is the unique global minimizer for all
`0<=g<=2.8934`.  Since `0<=A7<=A_full`, pointwise
`V7(p)>=V_full(p)`, so the same uniform-global statement transfers to rank 7
on that lower range.  No full-rank minimizer geometry is transferred.

For the upper side the package records a rational simplex witness with
denominator `10^9`.  Using the accepted outward lambda3--lambda6 intervals and
`mpmath.iv` on the exact rational probabilities gives

`R(p) in [3.7183448981203704..., 3.7183448981204046...]`

and at exact rational `g=3.71835`

`V_g(p) in [-2.3081604014e-6,-2.3081603858e-6]`.

Thus

`g_global in [2.8934,3.71835]`.

This does not identify the attaining orbit or prove a unique first transition.

## Stationary/minimum atlas and g=4

R7P-037 numerically saturates at three D12 orbit candidates for `g=3.7`, the
local-crossing gain, and `g=4`, and at fifteen candidates for `g=5`; the second
fixed-seed batch adds no orbit.  At `g=4` the observed stationary energies are
approximately `-0.1388476630`, `0`, and `0.0225772250`.  The first is therefore
the best **found** orbit, not a certified unique global orbit because the 7D
complement is not globally excluded.

The D12 action licenses orbit deduplication only.  Reflection-fixed C4,
two-harmonic, and fixed-amplitude phase reductions are invariant restricted
families, not global-minimizer forcing theorems.

## Declared gradient flow

For the explicitly chosen mathematical law `theta_dot=-grad Phi_4(theta)`,
the two numerical unstable branches from the observed index-one saddle go to
uniform and to the localized D12 orbit respectively.  This is a numerical
heteroclinic geometry for that Euclidean gradient flow only; no validated
trajectory tube, physical clock, mobility, nucleation rate, or basin weight is
claimed.

## Final M status

The rigorous global result is the gain bracket above.  Global orbit uniqueness
at `g=4`, stationary exhaustion, and validated dynamical connections remain
unresolved and are recorded as such rather than inferred from numerical
saturation.
