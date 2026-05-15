# P1751 / S701 — STRICT NONPROXY H1 BOUNDARY-CONTROL CONTRACT UPGRADE (PL)

Status: `P1751_EXECUTED_STRICT_NONPROXY_H1_BOUNDARY_CONTROL_CONTRACT_UPGRADE_NO_FALSE_PASS`

## Cel

Po `P1750` zaktualizować kontrakt nonproxy (`P1732`) tak, aby test H1 miał
jawnie **dwa poziomy raportowania**:

1. strict-local,
2. weak-form z klauzulą brzegową.

## Co zmieniamy

- Do kontraktu H1 dodajemy `boundary_control_contract` jako wymagany element
  przed decyzją PASS/OBSTRUCTION.
- Dodajemy schemat werdyktów:
  - `strict_local`: `PASS_ZERO` albo `OBSTRUCTION`,
  - `weak_form_with_boundary`: `PASS_WEAK_FORM_WITH_BOUNDARY_CLAUSE`
    albo `OBSTRUCTION`.
- Dodajemy zasadę promocji:
  weak-form pass **nie** promuje automatycznie strict-core closure.

## Znaczenie

To jest kolejny uczciwy krok w kierunku theorem-grade reverse chain:

`K_strict -> współczynniki -> pełny Lagrangian -> EOM -> H1(nonproxy)`

z formalną kontrolą wyrazów brzegowych, bez fałszywych passów.

## Następny uczciwy krok

Uruchomić pierwszy 4D nonproxy H1 `(A_μ,H)` zgodnie z tym kontraktem i zwrócić
równocześnie raport strict-local oraz weak-form.

## Plik artefaktu

- `generated/p1751_s701_strict_nonproxy_h1_boundary_control_contract_upgrade_checkpoint.json`
