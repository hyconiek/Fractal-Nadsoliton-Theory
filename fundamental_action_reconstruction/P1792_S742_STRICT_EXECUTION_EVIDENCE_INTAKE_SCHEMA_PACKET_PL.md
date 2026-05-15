# P1792 S742 Strict Execution Evidence Intake Schema Packet (PL)

## Cel

Dostarczyć jednolity schema-protokół przyjmowania wyników wykonania,
żeby `P1791` mogło działać deterministycznie na rzeczywistych artefaktach.

## Dlaczego to potrzebne

Mamy kontrakt (`P1790`) i protokół werdyktu (`P1791`),
ale nadal brakuje standardu wejściowego: jakie pola musi mieć każdy ledger,
aby uniknąć ręcznych interpretacji i niespójnych raportów.

## Intake schema (minimalny)

Każdy ledger dla `SV1..SV5` musi zawierać:

1. `ledger_id` (unikalny identyfikator),
2. `produced_by` (generator/runner id),
3. `background_family_id`,
4. `index_convention_id`,
5. `boundary_clause_id`,
6. `component_basis` (np. B1/B2/B3/C1/C2),
7. `result_kind` (`PASS_ZERO` lub `OPEN_OBSTRUCTION_WITH_TRACE`),
8. `residual_vector` (jawny, także dla PASS: wektor zerowy),
9. `obstruction_tags` (puste przy PASS),
10. `timestamp_utc`.

## Strict validation rules

- Brak któregokolwiek pola -> `INVALID_LEDGER_SCHEMA`.
- Niezgodność `background/index/boundary` między ledgerami -> `FREEZE_MISMATCH`.
- `PASS_ZERO` bez jawnego `residual_vector` == `INVALID_PASS_CLAIM`.

## Zastosowanie

Schema jest obowiązkowe dla:

- `L1_EA_componentwise_ledger`,
- `L2_EH_componentwise_ledger`,
- `L3_ELg_componentwise_ledger`,
- `L4_H1_residual_vector_ledger`,
- `L5_boundary_control_confirmation_ledger`.

## Następny uczciwy krok

Wygenerować checkpoint maszynowy, który formalizuje schema i mapowanie
`INVALID_*` / `FREEZE_MISMATCH` / `OPEN_OBSTRUCTION_WITH_TRACE` dla intake pipeline.

## Objaśnienie dla laika

To format „formularza dowodowego”: jeśli raport nie zawiera wymaganych pól,
nie może być użyty do zaliczenia testu, nawet gdy wygląda obiecująco.
