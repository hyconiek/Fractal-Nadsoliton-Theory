# P2716/S1666 inversion-odd pseudoscalar source acceptance audit

Status: `P2716_PSEUDOSCALAR_REPRESENTATION_ADMISSIBLE_BUT_NO_STRICT_SOURCE_VALUE`

## Representation-theoretic maps
- candidate_maps=4
- equivariant_maps=2

## Acceptance rows
- map={'-chi': '-omega', '+chi': '+omega'}: source_exported=False, fixes_now=False. The map would transfer a strict signed pseudoscalar value to the orientation torsor, but current artifacts export no nonzero inversion-odd source value.
- map={'-chi': '+omega', '+chi': '-omega'}: source_exported=False, fixes_now=False. The map would transfer a strict signed pseudoscalar value to the orientation torsor, but current artifacts export no nonzero inversion-odd source value.

## Decision
P2716 separates two issues.  Representation-theoretically, an inversion-odd pseudoscalar sign torsor can couple equivariantly to the +omega/-omega orientation torsor: two equivariant maps exist.  But current artifacts export no non-premise, nonzero signed pseudoscalar/chiral source value, so the maps have no strict input sign to transport and do not fix lambda or discharge QW-2191.

## Next honest step
A next admissible move must either construct/export a concrete strict inversion-odd pseudoscalar or chiral source value, then test one of the two equivariant maps as a bounded witness, or pivot to a different new typed object outside the closed lanes.  Without a signed source value, preserve the P2697-P2716 no-new-live-frontier certificate.
