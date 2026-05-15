# P1732 S682 Strict Nonproxy Export Requirements for H1 and Metric Residual Packet (PL)

Status: `P1732_EXECUTED_NONPROXY_EXPORT_REQUIREMENTS_CONTRACT`  
As of: `2026-05-15`

## Cel

Po `P1731` (obstruction) ustalić minimalny, jednoznaczny contract eksportów nonproxy,
który odblokowuje dwa krytyczne testy:

1. `H1`: `δE_A^μ/δH - δE_H/δA_μ`,
2. `metric`: `EL_g - E_{μν}` na bazie `B1/B2/B3/C1/C2`.

## Co wyeksportowano

1. Lista wymaganych eksportów nonproxy dla testu `H1`.
2. Lista wymaganych eksportów nonproxy dla residualu metrycznego.
3. Polityka decyzji dla obu testów: tylko `PASS_ZERO` albo `OBSTRUCTION`.

## Rygor

- strict-only,
- bez bridge do legacy,
- bez fałszywego pass,
- pełny lagranżian nieszkieletowy utrzymany jako anchor.

## Następny uczciwy krok (rekomendacja)

Wykonać eksport minimalnego pakietu nonproxy z `P1732`, potem natychmiast
uruchomić dwa testy (`H1`, `EL_g-E_{μν}`) i wydać wyniki bezpośrednio
w klasach `PASS_ZERO`/`OBSTRUCTION`.
