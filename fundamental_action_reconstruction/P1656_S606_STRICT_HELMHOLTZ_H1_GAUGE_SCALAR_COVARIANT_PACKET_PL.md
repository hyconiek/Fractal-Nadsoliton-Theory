# P1656 / S606 — Strict Helmholtz H1 gauge+scalar covariant-divergence witness

## Cel
Rozszerzyć lokalny witness H1 (P1655) na sektor gauge+scalar,
z jawnym operatorem dywergencji kowariantnej `D_μ`, aby przesunąć tor odwrotny
`EOM -> L_total` o kolejny krok theorem-level.

## Zakres strict-only
- `strict_only=true`, `legacy_bridge_used=false`,
- brak claimu globalnego domknięcia H1..H4,
- wynik wyłącznie `PARTIAL` (no-false-pass).

## Konstrukcja
- Bierzemy sektor:
  `L = -1/4 F_{μν}F^{μν} + (D_μ φ)^*(D^μ φ) - V(φ)`.
- EOM:
  - gauge: `D_μ F^{μν} = J^ν(φ,A)`,
  - scalar: `D_μD^μ φ + ∂V/∂φ^* = 0`.
- Lokalny test H1:
  symetria krzyżowej wariacji między kanałem `A_ν` i `φ`
  na poziomie części jawnie zależnej od sprzężenia minimalnego.

## Wyjście
- `generated/p1656_s606_strict_helmholtz_h1_gauge_scalar_covariant_summary.json`
