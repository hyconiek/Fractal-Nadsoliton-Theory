# P1489 — S4.39 QW-2191 Selector Source Candidate (PL)

Status: `P1489_EXECUTED_QW2191_SELECTOR_SOURCE_CANDIDATE_LOCAL_ONLY`
As of: `2026-05-13`

## Dlaczego ten krok?

Masz rację: sedno `QW-2191` to brak jawnego źródła wyboru kierunku.

Dlatego ten krok nie robi tylko testu „czy działa”, ale buduje **kandydata
źródła selektora** wprost na poziomie strict-only.

## Decyzja profesorska

Definiujemy lokalny kandydat źródła:

`S_src_v1 := Delta_SB / (|W_SM - W_GR| + eps)`

Interpretacja fizyczna:

- dodatni znak `S_src_v1` = preferencja kanału SM nad GR,
- moduł `|S_src_v1|` = siła orientacji,
- ograniczamy `|S_src_v1|` przez policy gate (`<= 1`) by uniknąć sztucznego
  „dopchnięcia” kierunku.

To jest uczciwa próba rozwiązania rdzenia problemu QW-2191.
