# P1985 / S935 — Strict ADM/Bianchi-I Non-GB Curvature-Squared Lapse Obstruction Witness

## Claim (licensed scope)
Using strict-only coefficients `a_R2, a_Ric2, a_Riem2` from the B1 backend,
the weighted non-GB ADM/Bianchi-I lapse operator

`EL_nonGB = a_R2*EL_R2 + a_Ric2*EL_Ric2 + a_Riem2*EL_Riem2`

is **not** identically zero in the anisotropic sector. Thus, after the P1984
Gauss-Bonnet cancellation witness, a genuine strict non-GB residual remains.

## Method
1. Import strict coefficient values from P1853 evaluated layer.
2. Import symbolic lapse Euler-operator differences from P1981/P1982/P1983.
3. Form weighted sum and simplify exactly in SymPy.
4. Verify anisotropic non-zeroness plus isotropic-limit zero and numeric replay.

## Verdict
- `EL_nonGB != 0` in anisotropic Bianchi-I witness lanes.
- Isotropic limit (`sigma_i = dsigma_i = d2sigma_i = 0`) returns zero.
- High-derivative terms persist in the anisotropic weighted residual.

## Not licensed
This packet does not prove ToE closure, background-independence closure,
Cutkosky closure, PO2/PO3 closure, or QW-2191 selector discharge.
