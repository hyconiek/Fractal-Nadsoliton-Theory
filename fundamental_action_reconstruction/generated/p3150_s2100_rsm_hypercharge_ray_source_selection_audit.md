# P3150/S2100 R_SM hypercharge-ray source-selection audit

Status: `P3150_RSM_HYPERCHARGE_RAY_SOURCE_SELECTION_CONDITIONAL_RAY_NO_STRICT_SOURCE`

## Constructed object
- `Y_SM^ray hypercharge consistency ray`
- Classification: `conditional_source_selection_witness_for_hypercharge_ratios_only`
- Scope: `derives the P3148 hypercharge ratios from Yukawa plus anomaly equations, not the absolute normalization or representation content`

## Finite theorem
`P3150_T1_hypercharge_ray_derivation_and_source_obstruction`: The one-family Yukawa-invariance plus anomaly equations form a 5 x 6 linear system of rank 5, so they select a one-dimensional hypercharge ray.  Normalizing the Higgs charge to 1/2 gives exactly the P3148 assignments q=1/6, u=-2/3, d=1/3, l=-1/2, e=1, h=1/2; the SU3^2 U1 and U1^3 anomaly checks vanish on the same ray.  This is a genuine conditional source-selection witness for hypercharge ratios, but it is not a strict nadsoliton source: it assumes the SM-like field content/Yukawa obligations, leaves absolute normalization free until a unit/charge convention is imposed, and does not select the representation content itself.

## Finite counts
- `linear_equation_rows`: `5`
- `unknowns`: `6`
- `matrix_rank`: `5`
- `nullity`: `1`
- `ray_witnesses_matching_P3148`: `1`
- `redundant_anomaly_checks_vanishing`: `1`
- `candidate_source_rows`: `4`
- `strict_accepted_source_rows`: `0`

## Normalized ray
- `q`: `1/6`
- `u`: `-2/3`
- `d`: `1/3`
- `l`: `-1/2`
- `e`: `1`
- `h`: `1/2`

## Candidate source rows
- `Y_SM^ray consistency constraints`: ray `True`, normalization `False`, representation content `False`, strict source `False`; conditional ray witness only: uses SM-like field/Yukawa/anomaly obligations
- `P3146 unit/action axioms`: ray `False`, normalization `False`, representation content `False`, strict source `False`; dimension scales do not select gauge charges
- `P1960/P1961 local Lie/BRST algebra`: ray `False`, normalization `False`, representation content `False`, strict source `False`; local algebra supplies consistency rules, not a unique matter registry
- `P3148 installed R_SM^ax registry`: ray `True`, normalization `True`, representation content `True`, strict source `False`; self-referential as a source: it is the installed target registry

## Decision
P3150 improves the P3148/P3149 axiom branch: the hypercharge ratios are not arbitrary once the one-family field pattern, Yukawa terms, and anomaly constraints are assumed; they are the unique ray of the finite system.

## Why this is not strict
The witness is conditional on an installed SM-like field pattern and Yukawa/anomaly obligations.  No current strict object selects that field pattern, the Higgs normalization, the charge unit, or the registry from nadsoliton data.

## Recommendation
Construct P3151 as a finite source-selection audit for the representation content itself: test whether any strict object selects the five chiral multiplet pattern plus Higgs dimensions (3,2), (bar3,1), (bar3,1), (1,2), (1,1), (1,2) without importing SM family data.  If no source exists, preserve R_SM^ax/Y_SM^ray as conditional phenomenology and pivot to unit-charge normalization or GR/EH nonproxy coupling.
