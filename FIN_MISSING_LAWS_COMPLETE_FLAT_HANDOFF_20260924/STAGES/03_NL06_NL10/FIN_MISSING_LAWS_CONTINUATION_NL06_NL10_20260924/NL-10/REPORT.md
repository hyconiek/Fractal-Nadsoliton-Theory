# NL-10 — nonlinear transport gate

Status: **PHASE_KINK_STOP__RADIAL_COHABITATION_DOMAIN_WALL_CANDIDATE_OPENED**

The phase-kink route is closed by NL-06, but the **full radial FIN potential**
contains a different structure that must not be confused with it.

At the certified first coexistence point

`g_eq = 3.7183448981203875`

the aligned C4 model has two inequivalent global minima:

- uniform: `s=0`;
- localized:
  `[1.81990358128 , 1.913989554669, 1.914569132547, 1.367203280196]`.

Their energies agree with zero to floating replay precision.

Along the deliberately simple reaction coordinate

`s(t)=t s_localized`, `0<=t<=1`,

the FIN potential is positive throughout the open interval and reaches

`Phi = 0.046635807553` at `t = 0.5127`.

Starting from that point and solving the **actual four-dimensional stationary
equations** finds the index-one saddle

`[0.940957067321, 1.001439402378, 0.962108853628, 0.686414983995]`

with energy `0.046555459529` and Hessian eigenvalues

`[-0.066994352962,  0.110005192342,  0.156502279211,  0.176911931409]`.

The straight path overestimates the true saddle barrier by only
`8.034802e-05` (0.173%).

This changes the matter-search logic.  The best current candidate for a
nonlinear wall is **not** a kink between translation phases.  It is a boundary
between the uniform and localized radial FIN phases at coexistence.

If one *adds* a scalar intercell stiffness for the reaction coordinate t,

`E[t]=int [kappa_t/2 (t_x)^2 + W(t)] dx`,
`W(t)=Phi(t s_localized,g_eq)`,

then W has two stable zeros and is positive between them.  The standard first
integral gives a static heteroclinic in this **declared reduced model**.  For
`kappa_t=1`, its dimensionless wall-tension integral is
`0.198877664625`.  Endpoint reaction-coordinate curvatures are
`1.084635186484` and `1.328689448718`.

But the straight radial line is not an exact invariant submanifold of the full
C4 FIN equations: the maximum orthogonal gradient along it is
`5.980770e-03`.  Therefore this is **not yet** a full four-component domain-
wall certificate.  More importantly, FIN has not yet sourced the intercell
stiffness tensor for the radial amplitudes or a reversible kinetic metric.

**Decision:** open a radial-domain-wall branch, but do not run moving-soliton
numerics yet.
