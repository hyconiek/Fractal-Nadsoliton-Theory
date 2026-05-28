# P2312/S1262 — strict lapse/shear source normalization attempt

- Status: `OPEN_OBSTRUCTION_NON_GB_LAPSE_SOURCE_NORMALIZATION_SIGN_INDEFINITE`
- Candidate source: `P1985_STRICT_NON_GB_LAPSE_RESIDUAL`
- Required lift: `0.0068`
- Self-energy lambda `||c||^2`: `0.0094608371836785`
- Direct non-GB lapse source sign-indefinite: `True`
- Same shear energy sign flip witnessed: `True`
- Strict provider-lift export: `False`
- G1/G3 update: `OPEN_UNCHANGED`
- Theorem fingerprint: `474ec20b8e2b197480f932cd5ec8c6b40edba1140aef17b96fa288b09fca8afb`

## Concrete audit result
The P1985 non-GB lapse residual is a real strict ADM/Bianchi-I object, but explicit witnesses with the same shear energy give opposite signs.  A direct fixed-sign normalization cannot be a monotone policy lift; absolute-value or positive-part fixes would add a convention/selector layer.

## Guardrail statement
P2312 does not close G1/G3, does not discharge QW-2191, does not add a selector premise, does not transfer legacy-kernel roles, and does not claim ToE closure.

## Next honest step
Either export a new strict dynamical constraint fixing the lapse/shear derivative orientation used by P1985, or stop the G1/G3 bridge route at the P2308-P2312 blocker stack and move to a genuinely new provider/source class.
