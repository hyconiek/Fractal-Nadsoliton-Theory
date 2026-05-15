# P1798 S748 Strict W1..W4 Full Export Completeness And BW Entry Guard Packet (PL)

Status: `P1798_EXECUTED_STRICT_W1_W4_FULL_EXPORT_COMPLETENESS_AND_BW_ENTRY_GUARD_PACKET_NO_FALSE_PASS`
As of: `2026-05-15`

## Cel

Zamknąć kolejny realny bottleneck wskazywany wcześniej (`W1..W4 not all delivered as FULL_EXPORT`):

1. zdefiniować twardą macierz kompletności `W1..W4`,
2. powiązać ją z wejściem do `G_BW`,
3. uniemożliwić wejście do BW na niepełnych eksporach nonproxy.

## Zakres strict-only

Brak legacy bridge, brak proxy degradacji, brak theorem-level claim.
To wyłącznie execution guard przed false-pass.

## Definicja kompletności W1..W4

Każdy element `Wk` musi mieć jednocześnie:

1. `artifact_id` (jawny export),
2. `classification.representation = NONPROXY`,
3. `classification.artifact_level = FULL_EXPORT`,
4. `freeze_id` zgodny z `freeze_id_common` runu,
5. `verification.check_command` + `verification.check_digest`.

## BW entry guard

`G_BW_ENTRY_ALLOWED = true` tylko gdy:

- `W1,W2,W3,W4` wszystkie mają status `FULL_EXPORT_VERIFIED`,
- brak freeze mismatch,
- brak wpisów `OPEN_OBSTRUCTION_WITH_TRACE` na poziomie W-matrix.

W przeciwnym wypadku:

- `G_BW_ENTRY_ALLOWED = false`,
- `G_BW = OPEN_OBSTRUCTION_WITH_TRACE`,
- `G_BRST/G_CUT = LOCKED`.

## Co jest dowiedzione

1. Bariery wejścia do BW są teraz jawnie sparametryzowane przez kompletność W-matrix.
2. Brak pełnych eksportów W1..W4 nie może być ukryty przez PASS z innego fragmentu pipeline.
3. Reguła jest zgodna z wcześniejszymi lockami (`P1767`, `P1795`, `P1797`).

## Co pozostaje OPEN

1. Rzeczywiste dostarczenie `FULL_EXPORT_VERIFIED` dla W1..W4 na wspólnym freeze.
2. Faktyczny `PASS_ZERO` dla BW po wejściu przez nowy guard.
3. Theorem-level BRST/Cutkosky closures.

## Produkt

- Packet PL.
- Checkpoint JSON + W-matrix template JSON do natychmiastowego wypełnienia.
