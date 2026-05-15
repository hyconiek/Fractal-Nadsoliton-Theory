# P1810 S760 Strict TG1 Lock Precedence Resolution Packet (PL)

Status: `P1810_EXECUTED_STRICT_TG1_LOCK_PRECEDENCE_RESOLUTION_PACKET_NO_FALSE_PASS`

## Cel

Usunąć proceduralną dwuznaczność `TG1_BW` między równoległymi audytami (`P1805`, `P1808`, `P1809`) przez jawny porządek priorytetu locków.

## Reguła priorytetu (od najwyższego)

1. `OPEN_LOCKED_BY_W1_CONSISTENCY_CONFLICT`
2. `OPEN_LOCKED_BY_UNIFIED_NONPROXY_RESIDUAL`
3. `PASS_ZERO`

Jeśli wiele źródeł podaje różne statusy, wybieramy status o wyższym priorytecie blokady.

## Co zostało dowiedzione

- Repo ma teraz jednoznaczną regułę scalania locków TG1 bez ręcznych decyzji.

## Co pozostaje OPEN

- Realne `PASS_ZERO` dla unified residual run-pack.

## Ryzyka false-pass

1. Użycie „łagodniejszego” locka przy istnieniu konfliktu W1.
2. Ręczne przestawienie TG1 bez przejścia przez precedence merger.

## Następny uczciwy krok

Po każdym nowym runie: uruchomić precedence merger i dopiero z tego artefaktu zasilać state-vector.
