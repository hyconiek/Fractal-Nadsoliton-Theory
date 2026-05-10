# P1179 Targeted Out-Of-Locality Uplift Scan Packet

Status: `P1179_EXECUTED_TARGETED_OUT_OF_LOCALITY_UPLIFT_SCAN_NO_FALSE_PASS`
As of: `2026-05-10`

## Goal

Wykonać następny uczciwy krok po `P1178`: znaleźć najlepszy lokalnie kandydat
pod kątem `robust_fraction` w out-of-locality scanie, bez zmiany twierdzeń
teoretycznych.

## Professor-level decision

Dodaję skan `p1179_targeted_out_of_locality_uplift_scan.py`, który uruchamia
`P1171` dla kandydatów sąsiednich (`p1170_neighbor_1..5`) i wybiera najlepszy
wynik do dalszych testów.

## Honest boundary

`P1179` pozostaje narzędziem selekcyjnym; brak claimu closure i brak `QW-2191`
discharge.
