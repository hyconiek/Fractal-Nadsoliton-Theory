# P1741 S691 Strict Full-Chain Forward/Reverse QG Dependency DAG Packet (PL)

Status: `P1741_EXECUTED_STRICT_DEPENDENCY_DAG_EXPORT`  
As of: `2026-05-15`

## Cel

Wyeksportować DAG zależności dla całego toru strict-only,
żeby było jawne, które kroki muszą poprzedzać aktualizacje bramek QG.

## Co wyeksportowano

1. Węzły i krawędzie DAG od `K_strict` do bramek QG.
2. Regułę: bez wyników `T1/T2` nie wolno uruchamiać `T3`.
3. Anchor pełnego lagranżianu nieszkieletowego.

## Rygor

- strict-only,
- bez bridge do legacy,
- bez fałszywego pass.

## Następny uczciwy krok (rekomendacja)

Wykonać `T1` i `T2` w kolejności z DAG, opublikować klasyfikacje,
a potem dopiero uruchomić `T3`.
