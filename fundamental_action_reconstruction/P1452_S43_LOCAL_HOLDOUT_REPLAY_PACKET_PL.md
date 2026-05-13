# P1452 S4.3 Local Holdout Replay Packet (PL)

Status: `P1452_S43_PROTOCOL_READY_STRICT_ONLY_NO_LEGACY_BRIDGE`
As of: `2026-05-13`

## Cel

Po `P1451 PASS_TARGETED_RERUN_LOCAL_ONLY` wykonać mały holdout lokalny (3 scenariusze), aby sprawdzić czy poprawka nie działa wyłącznie punktowo.

## Kontrakt S4.3

Dla każdego holdout case:

1. `gain >= min_gain`,
2. `replay_gap <= replay_tol`,
3. local-only scope i brak zamknięcia `QW-2191`.

Werdykty:

- `PASS_HOLDOUT_LOCAL_ONLY`
- `FAIL_HOLDOUT_MARGIN`
- `FAIL_HOLDOUT_REPLAY`
- `FAIL_HOLDOUT_SCOPE`

## Decyzja profesorska

To etap krytyczny: po lokalnym sukcesie na jednym przypadku sprawdzamy uogólnienie minimalne, ale nadal bez inflacji do global closure.

## Rekomendacja następnego uczciwego kroku

**Uruchomić P1452 i jeśli pojawi się FAIL — od razu eksportować holdout obstruction oraz wrócić do niecyklicznego projektu provider (S4.4).**

## Omówienie dla laika

To jak test leku na kilku nowych pacjentach po udanym przypadku pierwszym: jeden sukces to za mało, trzeba sprawdzić, czy efekt powtarza się choć w małej grupie.
