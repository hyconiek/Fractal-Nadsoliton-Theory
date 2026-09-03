# FIN ST2532–ST2621 — Configuration Generator and Dynamic-RG No-Go

## Executive conclusion

The strict FIN loop tensor and triangular-prism construction determine a
conditional equilibrium family and an exact stationary renormalization map.
They do not determine configuration-space kinetics.

For the 12-bit state space there are 4096 configurations and 24,576 undirected
single-spin-flip edges. Every reversible generator with stationary law
\(\pi\) has

\[
Q_{xy}=\frac{c_{xy}}{\pi_x},
\qquad c_{xy}=c_{yx}>0.
\]

Detailed balance fixes the forward/backward ratio but leaves one positive
activity \(c_{xy}\) per edge. After exact \(D_{12}\) quotienting, 1056 activity
orbits remain: 992 of size 24 and 64 of size 12.

## Kinetic counterexample

On one three-bit triangle at \(\theta=0.3\), Metropolis, heat-bath and symmetric
square-root rates have the same Gibbs state and satisfy detailed balance, but
their relaxation spectra differ.

| Rule | Spectral gap |
|---|---:|
| Metropolis | 1.31980017 |
| heat bath | 0.85213730 |
| symmetric | 1.78154388 |

The spectra are not related only by a common time rescaling. Thus stationary
state, locality and detailed balance do not determine dynamics.

## Maximum-caliber boundary

Maximum caliber does not remove the ambiguity. It requires:

- a reference path measure or prior generator,
- a mean activity/jump constraint,
- an Onsager mobility or edge conductance,
- and an overall time scale.

A uniform edge prior conditionally gives square-root rates, while other priors
give Metropolis- or heat-bath-like rules. Information optimization relocates
the kinetic assumptions rather than deriving them.

## Stationary RG versus dynamical RG

The symmetric triangular-prism Gibbs family has the exact stationary map

\[
\theta_{\rm eff}=\theta+rac12\log
\frac{e^\theta A(q)+e^{-\theta}B(q)}
     {e^{-\theta}A(q)+e^\theta B(q)}.
\]

However, projecting a six-bit prism dynamics onto the bottom triangle is not
strongly lumpable for generic \(q\neq0\). Two fine states with the same bottom
configuration but different top configurations give different bottom-flip
energy changes and hence different transition rates.

Therefore:

\[
\boxed{
\text{stationary RG closure does not imply dynamical RG closure}.}
\]

At \(q=0\), lumpability returns, but the refinement is dynamically trivial.
For \(q\neq0\), the coarse process generally acquires memory or must retain a
hidden top-layer state.

## Relation to dual FIN dynamics

The three generators act on different spaces:

- \(-iA\): unitary dynamics on 12 complex amplitudes,
- \(-A\): heat dynamics on 12 vertex values,
- \(Q\): stochastic dynamics on 4096 binary configurations.

No strict functor uniquely lifts \(A\) to \(Q\). Shared spectral organization
does not make the three operational channels identical. Configuration parity
is not reconstructible from vertex populations alone.

## Gate

| Requirement | Result |
|---|---|
| finite ternary Gibbs state | PASS |
| exact stationary prism RG | PASS |
| reversible generators exist | PASS |
| unique kinetic activity | FAIL |
| strict rate rule | FAIL |
| strict clock | FAIL |
| dynamic RG lumpability | FAIL |
| canonical \(A\to Q\) lift | FAIL |
| physical units | FAIL |
| OA evidence | FAIL |

Final score:

\[
\boxed{3/10\ \text{mathematical rows},\qquad0/10\ \text{physical rows}.}
\]

## Deepest surviving interpretation

FIN now supplies a strict interaction shape, a conditional non-Gaussian state
family, and an exact equilibrium prism-RG formula. It still does not supply the
kinetic activity, time scale, hidden-state handling, detailed-balance prior,
apparatus or physical calibration. Equilibrium information geometry is not yet
physical dynamics.

## Highest-value next theorem

Either derive a strict kinetic activity/clock source compatible with the prism
RG map, or prove an impossibility theorem for every local reversible
single-spin-flip generator sourced functorially from \(A\) and \(\tau\).

No selector, physical clock, continuum, laboratory evidence, Standard Model,
gravity, \(L_{\rm total}\), or Theory of Everything closure follows.
