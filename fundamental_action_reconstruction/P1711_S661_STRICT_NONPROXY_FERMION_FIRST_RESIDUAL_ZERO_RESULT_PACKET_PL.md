# P1711 S661 Strict Nonproxy Fermion First Residual-Zero Result Packet (PL)

Status: `P1711_EXECUTED_STRICT_FERMION_FIRST_RESIDUAL_ZERO_NO_FALSE_PASS`  
As of: `2026-05-15`

## Cel

Po `P1710` wykonać analogiczny pierwszy residual-zero test dla sektora fermionowego.

## Co wyeksportowano

1. Jawny `L_fermion_reduced`.
2. Równania `E_psi`, `E_psib`.
3. Residuale `EL - E` dla obu równań.
4. Wynik: `PASS_FERMION_ZERO`.

## Rygor

- strict-only,
- bez bridge do legacy,
- bez fałszywego pass,
- status globalny nadal `KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED`.

## Dla laika

To drugi zaliczony test jakości równań: tym razem dla fermionów. Oznacza to,
że nie tylko sektor gauge+Higgs, ale też sektor fermionowy jest lokalnie spójny
z lagranżianem. Nadal brakuje jednak sektora grawitacyjnego i dowodów globalnych.

## Następny uczciwy krok (rekomendacja)

Połączyć wyniki `P1710` i `P1711` w wspólny certyfikat częściowy, a potem
przejść do residual-zero dla sektora metrycznego i testów Bianchi/Ward.
