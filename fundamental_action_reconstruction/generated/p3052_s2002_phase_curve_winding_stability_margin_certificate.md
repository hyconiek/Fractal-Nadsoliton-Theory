# P3052/S2002 phase-curve winding stability-margin certificate

Status: `P3052_PHASE_CURVE_WINDING_STABILITY_MARGIN_CERTIFICATE_BOUNDED_NO_EXPORT`

## Finite certificate
- base integer winding: `1`
- clearance rows: `12`
- positive clearance rows: `12`
- minimum centroid-edge clearance: `0.04088112572302`
- perturbation rows: `40`
- winding-preserved perturbation rows: `40`
- Aut signed-winding rows: `4`
- Aut signed-winding verified rows: `4`
- source acceptance criteria: `4/8`
- satisfied proof obligations: `4/6`
- strict winding source theorem exported: `False`

## Decision
P3052 strengthens P3051 from a single winding witness to a finite stability-margin certificate: the centroid has positive clearance from every edge, all deterministic Fourier perturbation rows preserve winding, and Aut signed-winding behavior is explicit.  This is stronger evidence for a robust receiver-level orientation hint, but robustness is not a strict source theorem, selector/readout coupling, or L_total installation.

## Recommendation
Do not promote robust winding to selector closure.  The next proof-grade move must supply a strict source theorem explaining why this robust winding orientation is the nadsoliton selector sign, or pivot to an independent typed object outside sampled K/M geometry.
