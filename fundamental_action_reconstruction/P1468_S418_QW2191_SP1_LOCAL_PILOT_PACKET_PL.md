# P1468 — S4.18 QW-2191 SP1 Local Pilot A/B (PL)

Status: `P1468_EXECUTED_QW2191_SP1_LOCAL_PILOT_NO_GLOBAL_CLAIM`
As of: `2026-05-13`

## Cel

Wykonać minimalny test A/B dla kandydata `SP1_discrete_orientation_seed`:

- A: bez premisy SP1,
- B: z premisą SP1.

To tylko lokalny test feasibility; bez claimu strict-core closure.

## Kryteria

- wymagane: poprawa lokalnego wskaźnika selektora w B względem A,
- wymagane: jawne `strict_core_qw2191_closed=false`.

## Dyscyplina

- bez legacy bridge,
- bez global closure claim,
- status premisy pozostaje `NON_STRICT_UNLESS_PROVEN_INTERNAL`.
