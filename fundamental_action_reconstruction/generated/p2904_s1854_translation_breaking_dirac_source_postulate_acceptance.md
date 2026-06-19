# P2904/S1854 translation-breaking Dirac source postulate acceptance

Status: `P2904_TRANSLATION_BREAKING_DIRAC_SOURCE_POSTULATE_ACCEPTED_AS_CANDIDATE_NO_STRICT_PROVENANCE`

## Constructed source candidate
- `Xi_{0,+}`: signed Dirac source on `Z12`
- source values: `[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]`
- nonzero computed value: `1`
- translation stabilizer size: `1`
- translation orbit size: `12`

## Coupling
- selected basepoint: `0`
- selected polarity: `1`
- defect edge: `[0, 5]`
- coupling constructed: `True`

## Boundary
P2904 constructs the minimal translation-breaking object demanded by P2903: Xi_{0,+} has one nonzero signed value, trivial translation stabilizer, orbit size 12, and a unique coupling to A(0,+), D=(0,5), and the symbolic rho_9/5 template.  This passes the fixed-point obstruction as a candidate/source postulate.  It still does not prove strict nadsoliton provenance for why Xi_{0,+} rather than a translate or sign flip is exported, and it does not lift U_9_5 to a unit-bearing nonproxy L_total coupling.  Therefore closure remains blocked.

## Recommendation
The next proof-grade move should audit provenance of the new Xi_{0,+} source candidate: either derive its support and sign from a strict nadsoliton asymmetry/chiral/defect-generation theorem, or prove that no current artifact sources Xi rather than its 23 translated/sign-flipped alternatives.  Do not spend the next step on more translation-neutral selectors or more pointed templates.
