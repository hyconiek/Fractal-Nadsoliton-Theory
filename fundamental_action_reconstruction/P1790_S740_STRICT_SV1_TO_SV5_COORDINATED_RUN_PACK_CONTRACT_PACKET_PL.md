# P1790 S740 Strict SV1->SV5 Coordinated Run-Pack Contract Packet (PL)

## Cel

Zamknąć lukę wykonawczą wskazaną w `P1789`: uruchomić jeden skoordynowany
run-pack dla `SV1..SV5` (EA/EH/ELg + boundary + H1) w jednym standardzie wykonania.

## Kontrakt wykonania

- `strict_only`: true (zero legacy bridge).
- `nonproxy_only`: true.
- `single_background_family`: wymagane.
- `single_index_convention`: wymagane.
- `single_boundary_control_clause`: wymagane.
- `decision_policy`: tylko `PASS_ZERO` albo `OPEN_OBSTRUCTION_WITH_TRACE`.

## Run order (no-skip)

1. `SV1` — komponentowy `E_A^μ` export.
2. `SV2` — komponentowy `E_H` export.
3. `SV3` — komponentowy `EL_g` export + residual target binding.
4. `SV4` — boundary-term control check on same background/index freeze.
5. `SV5` — `H1`: `δE_A^μ/δH - δE_H/δA_μ`.

## PASS requirements

Aby wydać `PASS_ZERO` dla run-pack:

- jawny ledger komponentowy dla `E_A^μ`,
- jawny ledger komponentowy dla `E_H`,
- jawny ledger komponentowy dla `EL_g`,
- jawny ledger residual `H1` (zero vector),
- boundary-control potwierdzony na tym samym freeze.

Brak któregokolwiek elementu -> `OPEN_OBSTRUCTION_WITH_TRACE`.

## Co NIE wolno inferować

- `PASS_ZERO` run-packu **nie** oznacza theorem-level QG closure.
- Promotion `BW/BRST/Cutkosky` nadal podlega lockom z `P1788`.

## Następny uczciwy krok

Po wykonaniu run-packu zaktualizować state-vector:

- `SV1..SV5` na podstawie jawnych ledgerów,
- bez zmiany `SV6..SV8`, dopóki globalne witnessy nie są dostarczone.

## Objaśnienie dla laika

To jest dokładna instrukcja jednego „pełnego testu lokalnego” —
jeśli choć jeden wymagany wynik nie jest udokumentowany, test nie może być uznany za zaliczony.
