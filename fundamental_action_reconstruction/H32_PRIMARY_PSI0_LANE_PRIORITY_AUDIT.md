# H32 Primary psi0 Lane Priority Audit

Status: `PASS_PARTIAL_PRIMARY_LANE_PRIORITY_SET`
Date: `2026-03-06`

## Goal

Ustawic jawny porzadek badan po `H31` i `V1..V7`:

- lane `psi0` pozostaje glownym kandydatem anchoru,
- lane `informational viscosity` pozostaje lane wtornym,
- bez zamykania tej slabszej drogi.

## Inputs

- `H30`: `psi0` jest deterministycznym kandydatem anchoru z kernel invariants.
- `H31`: istnieje formalny embedding `psi0 -> pair1`.
- `V1..V7`: `informational viscosity` pozostaje fizycznie sensowna, ale tylko jako lane wtorny typu `anchor-amplifying / response-splitting`.

## Decision

Najbardziej uczciwy priorytet roboczy na obecnym stanie jest taki:

- `primary lane = psi0`
- `secondary lane = informational viscosity`

To oznacza:

- nowe kroki domykajace powinny najpierw testowac, czy `psi0` da sie uczciwie zredukowac do selektora,
- lane `viscosity` ma byc utrzymany jako konkurencyjna hipoteza pomocnicza,
- ale nie powinien juz byc traktowany jako rownorzedny kandydat primary selector source.

## Frontier

`H32_B1 := psi0 is now the primary working anchor candidate and informational viscosity is now the retained secondary lane, but no strict-core selector closure follows from this prioritization alone`

## Hard limits

- no `theorem-level PASS`
- no `full-closure PASS`
- no claim that `psi0` is already a strict-core selector datum
- no claim that `informational viscosity` is closed or discarded
- no claim that `QW-2191` is discharged
