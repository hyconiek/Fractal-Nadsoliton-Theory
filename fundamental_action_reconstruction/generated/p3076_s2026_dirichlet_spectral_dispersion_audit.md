# P3076/S2026 Dirichlet spectral-dispersion audit

Status: `P3076_INTERNAL_DIRICHLET_DIFFUSIVE_SPECTRAL_BRANCH_WAVE_LIGHTLIKE_OBSTRUCTION`

## Finite certificate
- content grep lanes: `3`
- content grep hits: `32757`
- P3075 local Dirichlet accepted rows: `96`
- Z12 modes: `12`
- nonconstant modes: `11`
- fractional steps: `3`
- amplification rows: `36`
- nonexpansive amplification rows: `36`
- strictly contracting nonconstant amplification rows: `33`
- oscillatory phase rows: `0`
- accepted lightlike branch rows: `0`
- satisfied proof obligations: `4/6`

## Decision
P3076 diagonalizes the P3075 local Dirichlet/Laplacian source on Z12 and finds an internal diffusive smoothing branch: the constant mode is neutral and all 11 nonconstant modes contract under the three audited fractional steps.  The spectrum has a formal lambda_j ~ k^2 small-k proxy, but the exported dynamics is first-order dissipative gradient flow, not a second-order, unit-bearing, Lorentzian, wave/lightlike equation.

## Recommendation
Construct one bounded second-order lift obstruction table for the same Z12 Dirichlet source: add the minimal candidate phase-space variables (rho, pi), symplectic form, and Hamiltonian H = 1/2*pi^2 + E_D, then test exactly which premises are new and whether they are internally sourced or merely imported.  Keep it scoped as a missing-source audit; do not promote it to observed light, gauge photons, spacetime EOM, units, empirical physics, L_total, bridge/role-transfer, selector closure, or ToE.
