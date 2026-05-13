# P1449 S4 Extended Perturbation Obstruction-First Packet (PL)

Status: `P1449_S4_PROTOCOL_READY_STRICT_ONLY_NO_LEGACY_BRIDGE`
As of: `2026-05-13`

## Cel

Po S3 przechodzimy do S4:

```text
extended perturbation stress test (local-only) with obstruction-first policy
```

Zasada: pierwszy kontrprzykład => natychmiastowy eksport obstruction, bez „ratowania” wyniku narracją.

## Kontrakt S4

Dla każdego scenariusza rozszerzonej perturbacji:

1. `margin_after >= margin_before + min_gain`,
2. `replay_gap <= replay_tol`,
3. `scope_of_pass == LOCAL_ONLY_NON_GLOBAL_CLAIM`,
4. `strict_core_qw2191_closed == false`.

Werdykty:

- `PASS_LOCAL_STRESS_ONLY`
- `FAIL_STRESS_MARGIN` (+ obstruction export)
- `FAIL_STRESS_REPLAY` (+ obstruction export)
- `FAIL_STRESS_SCOPE` (+ obstruction export)

## Decyzja profesorska

To jest rygorystyczny krok: zamiast szukać „jak przepchnąć PASS”, szukamy najpierw gdzie teoria lokalnie pęka pod trudniejszą perturbacją.

## Rekomendacja następnego uczciwego kroku

**Uruchomić P1449 na prerejestrowanym zbiorze stresowym i przy pierwszym naruszeniu zatrzymać pipeline z obstruction artifact; bez podejmowania globalnych wniosków.**

## Omówienie dla laika

To jak test samochodu na wybojach: nie wystarczy, że jedzie po gładkiej drodze.
Jeśli zawieszenie pęka na pierwszej dziurze, to uczciwy raport to „nieprzygotowane”, a nie marketingowy sukces.
