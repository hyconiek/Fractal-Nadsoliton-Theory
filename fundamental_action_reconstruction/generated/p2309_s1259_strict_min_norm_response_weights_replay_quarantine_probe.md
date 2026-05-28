# P2309/S1259 — min-norm response weights replay quarantine

- Status: `OPEN_OBSTRUCTION_TARGET_CALIBRATED_WEIGHTS_REPLAY_QUARANTINED`
- Required lift: `0.0068`
- Induced lift from exported weights: `0.0068`
- Replay all rows meet target: `True`
- Replay worst margin: `0.0005576116041714485`
- Strict admissible: `False` because weights are calibrated to the target lift rather than derived from strict dynamics.
- G1/G3 update: `OPEN_UNCHANGED`
- Theorem fingerprint: `033d14d8cf0769e1e193d70e80db21139374d50a035fd22cd26064d494672faf`

## Guardrail statement
P2309 does not close G1/G3, does not discharge QW-2191, does not add a selector premise, does not transfer legacy-kernel roles, and does not claim ToE closure.

## Next honest step
Replace the target-calibrated constraint w·c=lambda_star with a strict variational source for lambda_star or w_i from ADM/Bianchi-I equations; otherwise keep P2308/P2309 as blockers for G1/G3.
