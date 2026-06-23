# P3050/S2000 area-sign to P3046 coupling-polarity audit

Status: `P3050_AREA_SIGN_TO_P3046_POLARITY_COUPLING_AUDIT_BOUNDED_NO_EXPORT`

## Finite certificate
- area-sign rows: `6`
- nonzero area-sign rows: `5`
- neutral area rows: `1`
- coupling rows: `11`
- nonzero candidate coupling rows: `10`
- Aut-equivariant nonzero coupling rows: `10`
- polarity-selected rows: `0`
- accepted orientation-coupling rows: `0`
- source acceptance criteria: `3/7`
- satisfied proof obligations: `3/6`
- P3046 coupling polarity selected: `False`

## Decision
P3050 constructs the exact A_s-to-P3046 coupling-polarity audit left by P3049.  Nonzero area signs have the right inversion-odd type and every nonzero row admits Aut-equivariant maps to the P3046 polarity torsor, but there are always two opposite maps.  The neutral lag-6 row supplies no sign.  Therefore no unique nonconventional orientation/coupling-polarity theorem is exported.

## Recommendation
Do not replay A_s sign-to-polarity maps as selector closure.  This lane is bounded unless a genuinely new strict orientation-law/source object selects one of the two maps; otherwise pivot to a different typed object or preserve the P3048-P3050 no-export certificate.
