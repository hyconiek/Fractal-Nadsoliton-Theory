# P1535 S485 Strict Selector Tie-Break Candidate Packet (No Legacy Bridge)

Status: `P1535_EXECUTED_STRICT_SELECTOR_TIE_BREAK_CANDIDATE_PROVISIONAL`
As of: `2026-05-14`

## Cel

Następny uczciwy krok po `P1534`:

- rozwiązać remis gałęzi wykryty przez skaner `LEM_4`,
- zaproponować strict-only regułę tie-break,
- nadal bez claimu global closure.

## Zakres

`S485` wprowadza kandydata reguły tie-break:

1. jeśli dwie gałęzie mają równy `branch_physics_score`,
2. wybieramy gałąź o większym `provenance_depth`,
3. jeśli nadal remis, wybór przechodzi do statusu `unresolved_tie`.

To jest kandydat theoremu/reguły, nie finalny dowód.

## Kontrakt wyjścia

- `tie_detected` (bool),
- `tie_break_status` in `{resolved_by_provenance_depth, unresolved_tie, no_tie}`,
- `selected_branch_after_tie_break`,
- `lem4_status_update` in `{open, partial}`,
- `qw2191_closed=false`.

## PASS/FAIL

PASS jeśli remis jest jawnie wykrywany i rozstrzygany lub jawnie oznaczany jako
nierozstrzygnięty.

FAIL jeśli pipeline ukrywa remis albo automatycznie domyka `QW-2191`.
