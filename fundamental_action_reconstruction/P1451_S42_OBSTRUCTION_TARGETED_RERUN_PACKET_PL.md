# P1451 S4.2 Obstruction-Targeted Rerun Packet (PL)

Status: `P1451_S42_PROTOCOL_READY_STRICT_ONLY_NO_LEGACY_BRIDGE`
As of: `2026-05-13`

## Cel

Wykonać uczciwy rerun po P1450:

```text
rerun only on obstructed scenario (s4_case_02) with proposal-applied delta
```

bez rozszerzania claimów i bez legacy bridge.

## Kontrakt S4.2

1. Bierzemy wyłącznie scenariusz z pierwszego obstruction (`s4_case_02`).
2. Stosujemy `delta_margin_boost` z `P1450`.
3. Sprawdzamy ten sam próg `min_gain` i `replay_tol` co w P1449.
4. Wynik pozostaje local-only (`LOCAL_ONLY_NON_GLOBAL_CLAIM`).

Werdykty:

- `PASS_TARGETED_RERUN_LOCAL_ONLY`
- `FAIL_TARGETED_RERUN_MARGIN`
- `FAIL_TARGETED_RERUN_REPLAY`

## Decyzja profesorska

To minimalna, niecykliczna walidacja hipotezy naprawczej: jedna poprawka -> jeden scenariusz -> jeden werdykt.

## Rekomendacja następnego uczciwego kroku

**Jeśli PASS: przejść do S4.3 (mały holdout lokalny 2-3 scenariusze). Jeśli FAIL: eksportować obstruction v2 i projektować nowy provider proposal bez global claims.**

## Omówienie dla laika

Naprawiliśmy konkretny błąd, więc teraz sprawdzamy dokładnie ten sam przypadek.
To jak wymiana jednej części w silniku i krótki test, czy zniknął dokładnie ten konkretny problem.
