# P1500 — S4.50 QW-2191 Selector-Source + F→(L_SM+L_GR) Witness (PL)

Status: `P1500_EXECUTED_QW2191_SELECTOR_SOURCE_AND_F_MAPPING_WITNESS_LOCAL_ONLY`
As of: `2026-05-13`

## Cel

Zrealizować dwa brakujące warunki globalizacji z `P1499`:

- `R2`: eksport strict internal selector source,
- `R5`: eksport jawnego witness mapowania `F(nadsoliton) => L_SM + L_GR`.

## Decyzja profesorska

Definiujemy obiekty:

1. `S_internal_v1` — strict selector source object,
2. `W_Fmap_v1` — witness rozdziału `F` na kanał SM/GR pod tym samym selektorem.

Warunek fizyczny:

- spójność orientacji selektora,
- dodatniość kanałów SM/GR,
- brak legacy bridge.
