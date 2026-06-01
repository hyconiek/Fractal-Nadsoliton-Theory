# P2363 S1313: legacy->strict bridge moment transport into L_total

Status: `OPEN_PROGRESS_BRIDGE_COMPLETED_MOMENT_TRANSPORT_TO_LTOTAL_NO_SELECTOR_OR_ROLE_TRANSFER`

## Result

P2363/S1313 certifies bridge-completed moment transport: applying the explicit APD completion factor to the legacy kernel gives the same `c0,c1,c2` moments and effective `L_total` couplings as the strict kernel.

## Main Equalities

- `K_legacy*Q_APD - K_strict = 0`.
- Completed moment residuals: `{'c0': '0', 'c1': '0', 'c2': '0'}`.
- Completed effective-coupling residuals: `{'m2_eff': '0', 'lam_eff': '0', 'y_eff': '0', 'g_eff': '0', 'xi_eff': '0'}`.

## Negative Control

- Raw legacy moment max mismatch: `1.6349146124693021`.
- Amplitude-normalized legacy moment max mismatch: `0.30348403424955017`.

## Hard Limits

- No strict source theorem for APD is claimed.
- No selector premise or QW-2191 discharge is claimed.
- No legacy physical-role transfer is claimed.
- No full 4D EOM theorem closure or ToE closure is claimed.
