# P3058/S2008 unique-polarity coupling constraint verifier

Status: `P3058_UNIQUE_POLARITY_COUPLING_CONSTRAINT_VERIFIER_BOUNDED_NO_EXPORT`

## Finite certificate
- content grep lanes: `3`
- content grep hits: `2877`
- current constraints: `5`
- constraint intersections: `31`
- unique-polarity intersections: `0`
- all-current accepted polarities: `['minus_polarity', 'plus_polarity']`
- polarity-odd current constraints: `0`
- satisfied proof obligations: `3/5`

## Decision
P3058 attacks exactly one P3057 atom, new_unique_polarity_coupling.  The finite verifier constructs the missing atom as a two-torsor coupling and enumerates all 31 nonempty intersections of five current compatibility constraints.  Every current constraint is polarity-even, every intersection leaves both plus_polarity and minus_polarity, and zero intersections select a unique polarity.  The needed object is therefore sharper: a strict polarity-odd source-law boundary condition coupled to G_selector, not another polarity-blind compatibility check.

## Recommendation
Do not repeat polarity-even compatibility constraints.  The next proof-grade move should construct one explicit strict polarity-odd source-law boundary condition for G_selector, with a signed value that selects plus_polarity or minus_polarity non-premise and then proves compatibility with the P3057 row-import and square obligations; otherwise pivot to a different named P3057 atom while preserving no selector/readout/L_total/bridge/ToE export.
