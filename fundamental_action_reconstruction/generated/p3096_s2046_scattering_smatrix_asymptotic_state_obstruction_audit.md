# P3096/S2046 scattering/S-matrix asymptotic-state obstruction audit

Status: `P3096_SCATTERING_SMATRIX_ASYMPTOTIC_STATE_OBSTRUCTION_BOUNDED_NO_GO`

## Finite certificate
- content grep lanes: `3`
- content grep hits: `32682`
- P3095 accepted non-imported dispersion/propagating observable sources: `0`
- finite channel rows: `12`
- channel rows with in-state: `0`
- channel rows with out-state: `0`
- channel rows with asymptotic region: `0`
- Born transition amplitude rows: `144`
- transition rows with unit-normalized amplitude: `0`
- transition rows with cross-section semantics: `0`
- S-matrix unitarity proxy rows: `3`
- S-proxy rows with exact unitarity: `0`
- S-proxy rows with physical scattering operator: `0`
- cross-section proxy rows: `12`
- cross-section rows with area unit: `0`
- cross-section rows with detector semantics: `0`
- scattering candidates: `5`
- candidate gate rows: `40`
- accepted non-imported scattering/S-matrix sources: `0`
- satisfied proof obligations: `4/5`

## Decision
P3096 constructs the requested scattering/S-matrix asymptotic-state obstruction audit.  The Z12 Laplacian supports finite Fourier channel labels, Born-like transition amplitudes for a finite source potential, formal S=I+i g T unitarity-defect proxies, and formal cross-section proxies.  These are real scattering-like witnesses, but no internal artifact exports in/out asymptotic states, a physical unitary S-matrix, unit-normalized scattering amplitudes, empirical detector semantics, spacetime asymptotics, or a non-imported physical scattering source.  Imported continuum scattering and empirical detector templates pass only as imported templates.  Therefore no non-imported scattering/S-matrix source is exported.

## Recommendation
Pivot to exactly one adjacent standard-physics interface atom outside selector replay: construct a bounded thermodynamic-radiation/blackbody-spectrum obstruction audit for the Z12 Dirichlet/Laplacian branch, testing whether finite spectral mode counts, formal propagation/scattering proxies, partition weights, and energy-flux rows supply a non-imported radiation spectrum, temperature/energy unit calibration, photon/light interpretation, and empirical intensity readout without importing continuum statistical field theory, apparatus units, observed light, L_total, bridge/role-transfer, or ToE.
