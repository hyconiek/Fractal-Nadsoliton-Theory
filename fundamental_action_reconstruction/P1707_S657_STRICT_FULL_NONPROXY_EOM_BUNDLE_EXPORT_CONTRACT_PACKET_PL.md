# P1707 S657 Strict Full Nonproxy EOM Bundle Export Contract Packet (PL)

Status: `P1707_EXECUTED_STRICT_FULL_NONPROXY_EOM_CONTRACT_NO_FALSE_PASS`  
As of: `2026-05-15`

## Cel

Po `P1706` wyeksportować ścisły kontrakt pełnego nonproxy bundle EOM dla całego
łańcucha strict:

`kernel strict -> współczynniki -> pełny lagranżian -> pełne równania ruchu`.

## Co wyeksportowano

1. Jawny zestaw pól dla sektorów metryka/gauge/higgs/fermion.
2. Docelowe kontrakty równań EOM dla każdego sektora.
3. Kontrakt wymaganych obiektów wariacyjnych i tożsamości zgodności.
4. Kontrakt bidirectional (forward/reverse) + lista otwartych theorem-level kroków.

## Rygor

- strict-only,
- bez bridge do legacy,
- bez fałszywego pass,
- status globalny pozostaje `KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED`.

## Dla laika

To jest „specyfikacja finalnych równań” teorii w wersji bez uproszczeń. Dzięki temu
następny krok może już dostarczyć konkretne równania wprost z teorii, zamiast
kolejnej warstwy planu.

## Następny uczciwy krok (rekomendacja)

Wypełnić kontrakt jawnie: wyeksportować obliczalne formuły `E_{μν}`, `E_fermion`,
`E_gauge`, `E_higgs` na jednej spójnej konwencji matematycznej.
