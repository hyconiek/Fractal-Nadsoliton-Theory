# P1654 / S604 — Strict bidirectional theorem requirement matrix with quantum-gravity closure gates

## Cel
Następny uczciwy krok po P1653: sformalizować **theorem-level matrix** dla pełnego toru
`K_strict -> współczynniki -> L_total -> EOM` i toru odwrotnego,
z jawnie wydzielonymi bramkami ToE dla problemu kwantowej grawitacji.

## Zakres strict-only
- `strict_only = true`
- `legacy_bridge_used = false`
- brak transferu roszczeń legacy na kernel strict.

## Bramy QG wymagane do uczciwego domknięcia ToE
1. renormalizacja (kontrola UV/counterterm flow),
2. unitarność (brak ghostów + zgodność optyczna amplitud),
3. background independence (obserwable i dynamika bez uprzywilejowanego tła).

## Reguła no-false-pass
- jeśli dowolna brama theorem-level jest `OPEN`, status strict-core closure pozostaje `OPEN`.

## Wyjście
- `generated/p1654_s604_strict_bidirectional_theorem_requirement_matrix_with_qg_summary.json`
