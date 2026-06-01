# P2372 S1322: bridge-kernel direct band-polarity audit

Status: `OPEN_PROGRESS_BRIDGE_KERNEL_DIRECT_POLARITY_AUDIT_NO_QW2191_DISCHARGE`

## Result

P2372/S1322 checks whether the bridge-completed strict kernel itself supplies the P2371 band-polarity inequality for the d5 selector. It does not.

## Certificate

- Strict K1: `0.4699856726450201`.
- Strict K5: `0.02413122336363006`.
- Strict K1/K5: `19.476247248756152`; required `< 1/3`.
- Strict direct maximizer pair distribution: `{'4,0': 12}`.
- Legacy normalized K5 is negative: `-0.24649432866906726`.

## Hard limits

- This is a direct-pair-weight audit, not a theorem that K(d) is the selector action.
- No `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, or ToE closure is claimed.
