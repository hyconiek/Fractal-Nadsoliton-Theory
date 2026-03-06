# H8 MINIMAL COMPONENT CARRIER CONSTRUCTION SPEC

Status: `PASS_PARTIAL_COMPONENT_CARRIER_CONSTRUCTION_SPEC_READY`
As of: `2026-03-06`

## Goal

Zapisac minimalny construction/export spec potrzebny do przejscia od samego ansatzu
`K_obs` do jawnego pair-level carrieru, na ktorym mozna policzyc blok `2x2`
`A_1` dla `pair1 = (c1,s1)`.

## Methodological note

Poniewaz bazowy kernel zostal zbudowany bez sprzezenia `obs`, caly lane `H`
pozostaje lane rozszerzenia hipotezy operatorowej. Ten krok nie reinterpretowuje
obecnego strict core jako juz zawierajacego `K_obs`; zapisuje tylko minimalna forme,
ktora nowe carryery musialyby przyjac, aby hipoteza byla testowalna.

## Minimal admissible construction routes

Exactly one of the following routes must be instantiated to move beyond `H7`.

### Route A: direct composite export

Provide directly on `V_1 = span{c1,s1}` an exported real symmetric matrix

`A_1 = [[a_1, b_1], [b_1, d_1]]`

with declared basis `{c1,s1}` and declared provenance from the hypothesized operator chain.

This is the shortest admissible route.

### Route B: factored carrier chain

Declare explicit finite carriers and component actions:

- `E_1 : V_1 -> L_1`
- `G_light : L_1 -> L_1`
- `R_mat : L_1 -> M_1`
- `O_obs : M_1 -> M_1`

with:

- explicit finite bases for `V_1`, `L_1`, `M_1`,
- explicit matrix representatives or explicit symbolic action rules,
- explicit pullback rule producing
  `A_1 = P_1 E_1^* G_light^* R_mat^* O_obs R_mat G_light E_1 P_1`.

This is the longer but fully factorized admissible route.

## Minimal required declarations

For either route, the repo must contain all of:

1. declared current pair `pair1 = (c1,s1)`,
2. declared carrier `V_1 = span{c1,s1}`,
3. declared chosen construction route (`A` or `B`),
4. explicit exported carrier object(s),
5. explicit basis labels on every involved carrier,
6. explicit provenance note linking the carrier(s) to the `H3` ansatz,
7. explicit statement that no selector value is smuggled into the construction.

## Immediate reduction of the problem

`H8` reduces the current blocker from

> no component-action carrier exists

into the sharper question:

> has the repo instantiated either Route A or Route B for `pair1`?

## Best current conclusion

The repository now has a packet-ready construction spec for how component carriers
would have to appear. What remains open is no longer the shape of an admissible
construction, but the absence of any instantiated construction route for `pair1`.

## Frontier after H8

- `H8_B1 := no explicit chosen construction route (direct composite export A_1 or finite factored carrier chain) has yet been instantiated for pair1`
- `H7_B1 := no explicit component-action carrier exists for E_1, G_light, R_mat, O_obs on pair1 or V_1, and no exported composite representative A_1 is present` is reduced to construction-route absence level
- `T12_B1 := strict-core typing judgment with totality and uniqueness remains undischarged`
- `T2_B1 := bridge theorem still specified but not discharged`
- `C32_B2 := raw cross-pair overlap route remains degenerate`

## Anti-overclaim boundary

This spec does **not** show that:
- the light-feedback hypothesis is true,
- the light-feedback hypothesis is false,
- Route A exists,
- Route B exists,
- `(a_1,b_1,d_1)` have been computed,
- `QW-2191` is discharged,
- theorem-level or full-closure PASS has been achieved.
