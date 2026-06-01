# P2391 S1341: selector-epoch rebased bridge gap matrix

Status: `PASS_ONE_SELECTOR_EPOCH_BIT_REBASED_BETA_TORS_TRANSPORT_AND_ROLE_TRANSFER_STILL_OPEN`

## Result

P2391/S1341 rebases the older bridge component-gap matrix after the P1343/P1348 selector epoch and P2390 qualification.
The selector/source vector changes by Hamming distance `1`: only the `topological_phase_bit_chi11` generic selector bit flips from old gap wording to selector-present wording.
The rebased GF(2) matrix rank is `2`.

## Chi11 row after rebase

- Generic strict selector exported: `True`.
- Explicit `beta_tors -> chi11` transport: `False`.
- Role transfer allowed now: `False`.

## Hard limits

- This removes stale generic-selector-gap wording only; P2392 retires `beta_tors -> chi11` as an active selector-search target rather than proving that transport.
- No completed legacy-to-strict bridge, no legacy role-transfer theorem, no `L_total` promotion, no cap-density source theorem, and no ToE closure is claimed.
