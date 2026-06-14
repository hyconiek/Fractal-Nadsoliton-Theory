# P2725/S1675 chiral-bispectrum translation-flow signed-velocity no-go

Status: `P2725_TRANSLATION_FLOW_SIGNED_VELOCITY_NO_GO_NO_STRICT_DYNAMIC_CHIRAL_SOURCE`

## Finite dynamic test
- checked_flow_rows=264
- velocity_count=11
- nonzero_signed_velocity_count=0
Every finite translation-flow difference of the P2718 Im(B_{1,5}) marker is zero for both orientations and all 11 nonzero Z12 velocities.  A dynamic source built only from translation of this marker has no nonzero signed value to couple to the P2721 polarity pair.

## Acceptance
- accepted_as_strict_dynamic_chiral_source=False
- missing=['nonzero_signed_value_exported', 'coupled_to_p2721_polarity_pair', 'selects_exactly_one_polarity']

## Recommendation
Do not reuse translation flow of Im(B_{1,5}) as the missing source.  The next admissible proof-grade move must either introduce a new non-translation dynamic/chiral observable with a computable nonzero signed value and an explicit P2721 polarity-coupling theorem, or preserve the P2697-P2725 no-new-live-frontier certificate.
