# P1724 S674 Strict Componentwise Metric Residual Obstruction Map Packet (PL)

Status: `P1724_EXECUTED_STRICT_OBSTRUCTION_MAP_NO_FALSE_PASS`  
As of: `2026-05-15`

## Cel

Po `P1723` przejść od samego stubu residualu do mapy obstrukcji: które grupy
składników dokładnie blokują redukcję do zera.

## Co wyeksportowano

1. Grupowanie obstrukcji według typu składników tensorowych.
2. Linki pochodzenia składników do sektorów pełnego lagranżianu.
3. Priorytetowy następny krok redukcji bazy.

## Rygor

- strict-only,
- bez bridge do legacy,
- bez fałszywego pass,
- status globalny nadal `KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED`.

## Dla laika

To „mapa problemu”: wiemy dokładnie, które fragmenty rachunku grawitacyjnego są
niedopasowane i skąd biorą się w teorii. Dzięki temu kolejny krok jest dużo bardziej
precyzyjny.

## Następny uczciwy krok (rekomendacja)

Zredukować grupę pochodnych kowariantnych krzywizny do wspólnej bazy i ponownie
sprawdzić residual `ELg-Emunu`.
