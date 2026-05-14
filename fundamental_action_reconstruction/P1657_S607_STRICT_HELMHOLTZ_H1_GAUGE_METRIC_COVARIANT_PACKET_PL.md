# P1657 / S607 — Strict Helmholtz H1 gauge-metric covariant witness

## Cel
Kolejny uczciwy krok po P1656: dołączyć kanał metryczny i wyeksportować
lokalny witness H1 dla sprzężenia gauge-metric w torze odwrotnym `EOM -> L_total`.

## Zakres
- strict-only, bez bridge do legacy,
- wynik `PARTIAL` (bez fałszywego domknięcia),
- jawne utrzymanie otwartych bram QG: renormalizacja, unitarność, background independence.

## Konstrukcja
- Użyty sektor:
  `L = sqrt(-g)*[-1/4 g^{μα}g^{νβ}F_{μν}F_{αβ} + (M_Pl^2/2)R - Λ]`.
- EOM:
  - gauge: `∇_μ F^{μν} = J^ν`,
  - metric: `M_Pl^2 G_{μν} + Λg_{μν} = T_{μν}^{gauge} + ...`.
- Test H1 lokalny:
  zgodność krzyżowej wariacji kanałów `A_ν` i `g_{μν}`
  dla części Maxwell-Einstein na poziomie operatorowym.

## Wyjście
- `generated/p1657_s607_strict_helmholtz_h1_gauge_metric_covariant_summary.json`
