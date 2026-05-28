# P2313/S1263 — provider→margin route stop and new source-class contract

- Status: `OPEN_OBSTRUCTION_CURRENT_PROVIDER_MARGIN_ROUTE_STOPPED_NEW_SOURCE_CLASS_REQUIRED`
- All current route blockers active: `True`
- Current route stopped: `True` for the current P2300/P2281/P1985/P2310 interface class only.
- New source class required: `True`
- G1/G3 update: `OPEN_UNCHANGED`
- Theorem fingerprint: `29585181882ae1c91f376dce2d1aba3d7f2aa2cf4acd9ff83d76c0c1ae01f96b`

## Concrete audit result
P2308-P2312 now block the current route: no current response functional, target-calibrated weights are quarantined, self-energy lacks the theorem bridge, and the concrete non-GB lapse source is sign-indefinite.  Repeating this route without a new source class would be cyclic.

## Admission contract
Any P2314+ source must satisfy NSC1-NSC6: typed domain, strict internal origin, signed monotone response, selector discipline, replay/gate rows, and all guardrails.

## Guardrail statement
P2313 does not close G1/G3, does not discharge QW-2191, does not add a selector premise, does not transfer legacy-kernel roles, and does not claim ToE closure.

## Next honest step
Start a genuinely new source-class construction attempt satisfying NSC1-NSC6, or explicitly mark any Shannon/selector branch as non-strict until it exports an internal selector/orientation source and survives P2302/P2282 replay.
