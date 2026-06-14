# P2711/S1661 inversion-odd character source-law sign-coupling audit

Status: `P2711_SOURCE_LAW_SIGN_COUPLING_DEGENERACY_NO_STRICT_SIGN_SOURCE`

## Finite candidate count
- anti_inversion_characters=2
- source_law_candidates=4
- sign_degenerate_pairs=2

## Degeneracy rows
- `chi_5_+1_7_-1_11_-1`: lambda -> -lambda exchanges +omega and -omega; strict_sign_source_exported=False
- `chi_5_-1_7_+1_11_-1`: lambda -> -lambda exchanges +omega and -omega; strict_sign_source_exported=False

## Decision
P2711 enumerates the four finite source-law candidates obtained from two inversion-odd characters times two coupling signs.  Each character has a lambda -> -lambda degeneracy that exchanges +omega and -omega.  Current artifacts export no non-premise law fixing lambda, so the candidate family remains a premise-sign pair rather than a strict selector source.

## Next honest step
A further admissible move must introduce an actually new strict mechanism fixing the coupling sign lambda, or pivot outside the selector/sign lane to a new typed object.  Without such a mechanism, preserve the P2697-P2711 no-current-unlock certificate and do not promote QW-2191, pair12 strict-core, L_total, role transfer, bridge closure, or ToE.
