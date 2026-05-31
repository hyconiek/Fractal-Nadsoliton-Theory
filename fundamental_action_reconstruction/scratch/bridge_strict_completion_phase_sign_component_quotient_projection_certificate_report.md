# Phase-sign component-quotient projection certificate

Status: `component-quotient-projection-section-and-boundary-restriction-certified-over-gf2`

## Result

This certificate exports the representative projection Q and boundary/interior
selector matrices and checks Q*S=I_5, S*Q projector idempotence,
G*B_path*S=B_quotient, and H*B_path*S=0 over GF(2).

## Summary

- Q projects node bits to component bits: `True`
- Q*S is identity: `True`
- S*Q fixes audited node bits: `True`
- Boundary restriction equals quotient coboundary: `True`
- Interior restriction vanishes: `True`
- Rank Q: `5`
- Rank S*Q: `5`

## Hard limits

- K_strict_gate remains the current live/full operational kernel.
- No unqualified identity K_legacy_ont == K_strict_gate is claimed.
- No proof derives A(d), P(d), D(d), omega/phi, beta/eta, or the transport cocycle from strict nadsoliton dynamics.
- No beta_tors -> chi_11 theorem is claimed.
- No legacy physical-role transfer to K_strict_gate is claimed.
- No QW-2191 selector discharge is claimed.
- No ToE closure is claimed.
