# P1715 S665 Strict Metric Curvature-Term Expansion Scaffold Packet (PL)

Status: `P1715_EXECUTED_STRICT_METRIC_CURVATURE_EXPANSION_NO_FALSE_PASS`  
As of: `2026-05-15`

## Cel

Po `P1714` rozpisać jawnie składniki krzywiznowe `H^(R2)`, `H^(Ric2)`, `H^(Riem2)`
w równaniu metrycznym.

## Co wyeksportowano

1. Jawne formuły scaffold dla trzech tensorów krzywiznowych.
2. Złożone równanie metryczne z tymi składnikami i `T_matter`.
3. Lista dalszych kroków do testu residual-zero sektora metrycznego.

## Rygor

- strict-only,
- bez bridge do legacy,
- bez fałszywego pass,
- status globalny nadal `KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED`.

## Dla laika

To etap „rozpisania trudnych członów grawitacji”. Dzięki temu możemy w kolejnym
kroku naprawdę sprawdzić, czy równanie grawitacyjne pochodzi z teorii bez ukrytych
skrótów.

## Następny uczciwy krok (rekomendacja)

W tej samej konwencji indeksowej policzyć jawny test `EL_g - E_munu = 0`
z użyciem tych rozpisanych członów.
