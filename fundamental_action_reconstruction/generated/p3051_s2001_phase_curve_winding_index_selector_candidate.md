# P3051/S2001 phase-curve winding-index selector candidate

Status: `P3051_PHASE_CURVE_WINDING_INDEX_SELECTOR_CANDIDATE_BOUNDED_NO_EXPORT`

## Finite certificate
- base integer winding: `1`
- base nonzero integer winding: `True`
- Aut/translation rows: `48`
- translation-stable rows: `48`
- orientation-reversing rows: `24`
- strict source exported rows: `0`
- source acceptance criteria: `4/8`
- satisfied proof obligations: `3/6`
- strict winding selector source exported: `False`

## Decision
P3051 constructs a genuinely different finite object: the global winding/turning index of the (K,M) phase curve.  The base curve has winding +1, all translations preserve the signed value, and Aut inversion units reverse it.  This is a real topological orientation hint, but it remains receiver-level because K/M provenance, nonconventional orientation law, P3046 coupling, selector/readout installation, and unit-bearing action/EOM are absent.

## Recommendation
Do not promote winding +1 as selector closure.  A next proof-grade move needs a strict source theorem for the phase-curve orientation/winding sign or an independent typed object outside sampled K/M receiver geometry; otherwise preserve the P3048-P3051 bounded no-export certificate.
