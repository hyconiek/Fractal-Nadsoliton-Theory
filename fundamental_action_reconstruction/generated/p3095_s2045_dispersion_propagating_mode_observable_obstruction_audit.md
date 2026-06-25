# P3095/S2045 dispersion/propagating-mode observable obstruction audit

Status: `P3095_DISPERSION_PROPAGATING_MODE_OBSERVABLE_OBSTRUCTION_BOUNDED_NO_GO`

## Finite certificate
- content grep lanes: `3`
- content grep hits: `32465`
- P3094 accepted non-imported stress-energy/metric-response sources: `0`
- dispersion/group-velocity rows: `12`
- dispersion rows with spacetime speed: `0`
- dispersion rows with observed-light semantics: `0`
- mode-packet evolution rows: `12`
- packet rows with physical time unit: `0`
- packet rows with detector observable: `0`
- Green-pole catalog rows: `48`
- Green-pole rows with asymptotic state: `0`
- Green-pole rows with scattering readout: `0`
- energy-flux proxy rows: `12`
- energy-flux rows with physical flux unit: `0`
- energy-flux rows with radiation interface: `0`
- propagation candidates: `5`
- candidate gate rows: `40`
- accepted non-imported dispersion/propagating observable sources: `0`
- satisfied proof obligations: `4/5`

## Decision
P3095 constructs the requested dispersion/propagating-mode observable obstruction audit.  The Z12 Laplacian supports a finite dispersion curve, group-velocity proxies, formal dimensionless mode-packet evolution, mass-regularized Green-pole catalogs, and modal energy-flux proxies.  These are real propagation-like witnesses, but no internal artifact exports a physical propagating field mode, spacetime speed/light cone, detector-independent observable, observed-light/radiation interface, physical time/apparatus units, or a non-imported physical propagation source.  Imported continuum wave and observed-radiation templates pass only as imported templates.  Therefore no non-imported dispersion/propagating observable source is exported.

## Recommendation
Pivot to exactly one adjacent standard-physics interface atom outside selector replay: construct a bounded scattering/S-matrix asymptotic-state obstruction audit for the Z12 Dirichlet/Laplacian branch, testing whether Green-pole catalogs, finite mode packets, source/response kernels, and propagation proxies supply non-imported in/out states, an S-matrix or cross-section, unit-normalized scattering amplitudes, and empirical detector semantics without importing continuum scattering theory, spacetime asymptotics, apparatus units, L_total, bridge/role-transfer, or ToE.
