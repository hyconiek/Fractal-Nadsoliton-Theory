# P2872/S1822 cyclic-predecessor endpoint orientation no-source audit

Status: `P2872_CYCLIC_PREDECESSOR_ENDPOINT_ORIENTATION_NO_SOURCE_AUDIT_NO_CLOSURE`

## Cyclic predecessor audit
- candidate class: `ambient Z12 cyclic-order predecessor localizer d=-1 mod 12, with reflection-symmetry acceptance test`
- predecessor endpoint: `11`
- exact predecessor target representation: `True`
- reflection-invariant adjacent predicates: `[{'selected': [], 'reflected_selected': [], 'reflection_invariant': True, 'selects_predecessor_11': False}, {'selected': [1, 11], 'reflected_selected': [1, 11], 'reflection_invariant': True, 'selects_predecessor_11': False}]`

## Boundary
P2872 shows that the cyclic predecessor convention d=-1 mod 12 exactly represents the endpoint target, but it is orientation-sensitive: reflection swaps d=1 and d=11, and current artifacts do not export a strict orientation/boundary-arrow source law or unit-bearing coupling theorem.

## Recommendation
Do not replay cyclic predecessor, d=-1, or clockwise/counterclockwise endpoint conventions as sourcehood.  A next proof-grade move must export a strict non-premise orientation/boundary-arrow law that fixes predecessor over successor and supplies the unit-bearing 9/5 coefficient/coupling theorem, or pivot to a different new typed object; otherwise preserve no-new-live-frontier.
