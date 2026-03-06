# H11 MINIMAL ROUTE A PROVENANCE SPEC

Status: `PASS_PARTIAL_ROUTE_A_PROVENANCE_SPEC_READY`
As of: `2026-03-06`

## Goal

Zapisac minimalny provenance spec potrzebny do tego, aby `A_1_cand` przestalo byc
samym carrier placeholder i moglo byc traktowane jako provenance-valid export na lane
rozszerzenia hipotezy operatorowej.

## Methodological note

Poniewaz bazowy kernel nie zawieral sprzezenia `obs`, provenance dla `A_1` musi jawnie
stwierdzac, ze chodzi o lane rozszerzenia operatorowego. Kazdy zapis, ktory ukrywa ten fakt,
bylby metodologicznie niepoprawny i liczylby sie jako selector smuggling lub retroaktywna
reinterpretacja strict core.

## Minimal required provenance fields

A provenance-valid `Route A` export for `pair1` must contain all of the following:

1. `object_id = A_1`
2. `basis = {c1,s1}`
3. `carrier = V_1 = span{c1,s1}`
4. `lane = hypothesis_extension_only`
5. `base_kernel_contains_obs = false`
6. `construction_mode = direct_composite_export`
7. `operator_origin = one_of{ exported_composite_A_1 , pullback_from_E_1_G_light_R_mat_O_obs }`
8. `selector_smuggling = false`
9. `strict_core_reinterpretation = false`
10. `coefficient_status = one_of{ unresolved , partially_exported , fully_exported }`
11. `provenance_gap_statement` if coefficients are not fully exported
12. `hard_limit_theorem_level_pass = false`
13. `hard_limit_full_closure_pass = false`

## Minimal template

A packet-ready minimal template is:

```json
{
  "object_id": "A_1",
  "basis": ["c1", "s1"],
  "carrier": "V_1 = span{c1,s1}",
  "lane": "hypothesis_extension_only",
  "base_kernel_contains_obs": false,
  "construction_mode": "direct_composite_export",
  "operator_origin": "UNRESOLVED",
  "selector_smuggling": false,
  "strict_core_reinterpretation": false,
  "coefficient_status": "unresolved",
  "provenance_gap_statement": "UNRESOLVED_FROM_OPERATOR_EXPORT_CHAIN",
  "hard_limit_theorem_level_pass": false,
  "hard_limit_full_closure_pass": false
}
```

## Immediate reduction of the problem

`H11` reduces the current question from:

> what counts as provenance-valid Route A?

into:

> does the repo contain any populated provenance instance matching this minimal spec?

## Best current conclusion

The repository now has a packet-ready minimal provenance spec for `Route A`. What remains open is not the shape of admissible provenance, but the absence of a populated provenance instance satisfying it.

## Frontier after H11

- `H11_B1 := no populated provenance-valid Route A instance exists yet for pair1 even though the minimal provenance spec is now packet-ready`
- `H10_B1 := a minimal Route A candidate instance exists, but no provenance-valid exported A_1 derived from the operator chain exists yet` is reduced to provenance-instance absence under an explicit spec
- `T12_B1 := strict-core typing judgment with totality and uniqueness remains undischarged`
- `T2_B1 := bridge theorem still specified but not discharged`
- `C32_B2 := raw cross-pair overlap route remains degenerate`

## Anti-overclaim boundary

This spec does **not** show that:
- a provenance-valid Route A export exists,
- `A_1` is derived from the operator chain,
- coefficients have been exported,
- Route A is discharged,
- Route B is unnecessary,
- the light-feedback hypothesis is true,
- `QW-2191` is discharged,
- theorem-level or full-closure PASS has been achieved.
