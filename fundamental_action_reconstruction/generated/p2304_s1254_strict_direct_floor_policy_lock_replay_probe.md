# P2304/S1254 — strict direct-floor policy-lock replay

- Status: `OPEN_OBSTRUCTION_STRICT_DIRECT_FLOOR_POLICY_LOCK_REPLAY_REFUTED`
- P2302 required lift: `0.0068`
- P2303 strict direct floor: `0.005223026434758646`
- Direct-floor shortfall: `0.0015769735652413535`
- Direct-floor replay all rows meet target: `False`
- Direct-floor worst margin: `-0.05087095982439993`
- Strict closure status: `HELD_OPEN_DIRECT_FLOOR_ROUTE_REFUTED`
- Theorem fingerprint: `dd6eaefcc08fe28e6df36b11fea13467bdec346a2eea8670a6182cb872192c09`

## Guardrail statement
P2304 refutes the conservative direct-floor closure route. It does not close G1/G3, does not discharge QW-2191, does not add a selector premise, does not transfer legacy-kernel roles, and does not claim ToE closure.

## Next honest step
Attempt a strict aggregation/norm-to-margin theorem for the P2300 operator basis; if that cannot be proven without a selector premise, keep G1/G3 open and route the sufficient positive/norm bounds as non-strict controls only.
