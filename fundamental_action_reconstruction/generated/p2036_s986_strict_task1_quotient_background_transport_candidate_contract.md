# P2036 S986 Strict Task-1 Quotient Background Transport Candidate Contract

Status: `OPEN_CANDIDATE_CONTRACT_WITH_TRACE`

Result kind: `OPEN_NEW_PREMISE_BACKGROUND_TRANSPORT_CONTRACT_CANDIDATE_EXPORTED__NO_GLOBALIZATION`

## Professor Decision

`EXPORT_NONSTRICT_BACKGROUND_TRANSPORT_CONTRACT_CANDIDATE_KEEP_TASK1_LOCAL`

P2036 exports the contract shape needed after P2035.  It is explicitly a
new-premise candidate, not a strict transport theorem.

## Candidate Contract

Contract id: `Task1_B1_quotient_background_transport_candidate_v1`

Source object: `P2034 local scalar B1 quotient class [a]_B1`

Quotient basis:

`['R2_bar', 'Ric2_bar', 'Riem2_bar']`

Null direction:

`(1, -4, 1, -1)`

Zeroth-order seed:

`M_FRW_from_B1 = I3`, `M_BianchiI_from_B1 = I3` at `sigma2=0`

This identity seed is only a boundary seed.  It is not background
independence.

## Acceptance Gates

- `C1_QUOTIENT_SOURCE`: `SATISFIED_BY_P2034_LOCAL_SOURCE`
- `C2_BASIS_PRESERVING_MAP`: `OPEN_SYMBOLIC_PLACEHOLDER`
- `C3_FINITE_PART_SCHEME_TRANSPORT`: `OPEN_SYMBOLIC_PLACEHOLDER`
- `C4_GB_NULL_TRANSPORT`: `OPEN_TOPOLOGICAL_CLASSIFICATION_MISSING`
- `C5_ANISOTROPIC_WITNESS_DATA`: `OPEN_LOOP_DATA_MISSING`
- `C6_GLOBAL_COCYCLE`: `OPEN_COCYCLE_THEOREM_MISSING`

## Not Licensed

- background-global Task-1 renormalization
- FRW or BianchiI transported counterterm values
- finite-part scheme transport theorem
- GB null direction topological transport theorem
- global atlas cocycle theorem
- tensor-component B1 renormalization
- QW-2191 selector closure
- ToE closure

## Next Honest Step

Compute or import one concrete same-scheme finite-part map on the
`R2_bar/Ric2_bar/Riem2_bar` quotient basis, then test C2-C3.
