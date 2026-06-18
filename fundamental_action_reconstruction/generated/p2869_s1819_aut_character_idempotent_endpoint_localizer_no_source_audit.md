# P2869/S1819 Aut-character idempotent endpoint-localizer no-source audit

Status: `P2869_AUT_CHARACTER_IDEMPOTENT_ENDPOINT_LOCALIZER_NO_SOURCE_AUDIT_NO_CLOSURE`

## Exact idempotent representation
- candidate class: `Aut(Z12)-character Fourier idempotent on the unit orbit {1,5,7,11}, scaled by 9/5`
- endpoint weights: `{'1': {'fraction': '0/1', 'numerator': 0, 'denominator': 1, 'denominator_prime_factors': {}, 'z12_compatible_denominator': True}, '5': {'fraction': '0/1', 'numerator': 0, 'denominator': 1, 'denominator_prime_factors': {}, 'z12_compatible_denominator': True}, '7': {'fraction': '0/1', 'numerator': 0, 'denominator': 1, 'denominator_prime_factors': {}, 'z12_compatible_denominator': True}, '11': {'fraction': '1/1', 'numerator': 1, 'denominator': 1, 'denominator_prime_factors': {}, 'z12_compatible_denominator': True}}`
- scaled prime vector: `{'2': {'fraction': '0/1', 'numerator': 0, 'denominator': 1, 'denominator_prime_factors': {}, 'z12_compatible_denominator': True}, '3': {'fraction': '0/1', 'numerator': 0, 'denominator': 1, 'denominator_prime_factors': {}, 'z12_compatible_denominator': True}, '5': {'fraction': '0/1', 'numerator': 0, 'denominator': 1, 'denominator_prime_factors': {}, 'z12_compatible_denominator': True}, '7': {'fraction': '0/1', 'numerator': 0, 'denominator': 1, 'denominator_prime_factors': {}, 'z12_compatible_denominator': True}, '11': {'fraction': '9/5', 'numerator': 9, 'denominator': 5, 'denominator_prime_factors': {5: 1}, 'z12_compatible_denominator': False}}`
- coefficient denominator support: `[2, 5]`

## Boundary
P2869 shows that Aut-character Fourier idempotents can represent the endpoint exactly, but only by importing the endpoint/polarity projector and the scaled 9/20 character coefficients.  This is representability, not strict sourcehood or a unit-bearing coupling theorem.

## Recommendation
Do not replay Aut-character idempotent endpoint projectors as sourcehood.  A next proof-grade move must supply a strict non-premise law selecting the character projector/polarity and coupling it with a unit-bearing coefficient theorem, or pivot to a different genuinely new typed object; otherwise preserve no-new-live-frontier.
