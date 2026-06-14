# P2714/S1664 Z12 orientation-torsor global-section obstruction

Status: `P2714_ORIENTATION_TORSOR_GLOBAL_SECTION_OBSTRUCTION_NO_STRICT_LAMBDA_FIX`

## New typed candidate
- `orientation_torsor_of_P2708_boundary_cocycle_line` with fiber `['+omega', '-omega']`

## Aut action
- unit `1` maps `+omega` to `+omega`; reversing=False
- unit `5` maps `+omega` to `+omega`; reversing=False
- unit `7` maps `+omega` to `-omega`; reversing=True
- unit `11` maps `+omega` to `-omega`; reversing=True

## Global section candidates
- `-omega`: aut_compatible=False, failure_count=2
- `+omega`: aut_compatible=False, failure_count=2

## Decision
P2714 tests one new typed candidate admitted after P2713: the orientation torsor of the P2708 boundary-cocycle line.  The finite Aut(Z12) action has orientation-reversing units 7 and 11, so neither +omega nor -omega is an Aut-compatible global section.  The torsor is real, but it remains a two-point premise-sign torsor rather than a strict mechanism fixing lambda.

## Next honest step
Do not replay the orientation-torsor, character, lambda-coupling, damping-to-selector, older-release, or direct-route lanes.  A next admissible move must supply a new strict source that breaks the orientation torsor by an internal non-premise law, or a different genuinely new typed object outside the closed lanes.  Without that, preserve the P2697-P2714 no-new-live-frontier certificate.
