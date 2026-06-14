# P2710/S1660 finite Aut(Z12) anti-inversion orientation-character source test

Status: `P2710_AUT_Z12_ANTI_INVERSION_CHARACTER_EXISTS_BUT_NO_STRICT_SOURCE`

## Character table summary
- total_characters=4
- anti_inversion_characters=2

## Candidate source rows
- `chi_5_-1_7_+1_11_-1`: strict_source_available=False, exports_orientation_source=False. The character is a possible parity label, not an exported strict law choosing +omega over -omega. Current artifacts provide no non-premise source selecting this character and its sign.
- `chi_5_+1_7_-1_11_-1`: strict_source_available=False, exports_orientation_source=False. The character is a possible parity label, not an exported strict law choosing +omega over -omega. Current artifacts provide no non-premise source selecting this character and its sign.

## Decision
The finite character table of Aut(Z12) contains exactly two inversion-odd characters, so the mathematical parity labels exist.  However current strict artifacts do not export either character as a non-premise physical/source law, and an abstract character still does not choose +omega rather than -omega.  Thus P2710 does not discharge QW-2191 or unlock pair12 strict-core, L_total, role transfer, bridge closure, or ToE.

## Next honest step
The next admissible move must supply a genuinely new strict artifact that couples one inversion-odd character to the boundary-cocycle sign with a non-premise sign convention, or pivot to a different new typed object.  Without that new source, preserve the P2697-P2710 no-current-unlock certificate.
