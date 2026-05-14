# P1661 / S611 — Strict H2 boundary-flux lemma scaffold (gauge-scalar)

## Cel
Po P1660: wyeksportować theorem-level scaffold lematu
`candidate asymptotics => boundary flux vanishing` dla `J^μ` w bloku gauge-scalar.

## Zakres
- strict-only, bez legacy bridge,
- wynik `PARTIAL` (szkielet dowodu + brakujące lematy analityczne),
- final strict-core closure pozostaje `OPEN`.

## Konstrukcja
- Teza lematu (robocza):
  dla pól z klasy `A_μ=O(r^{-1-ε}), D_μ φ=O(r^{-2-ε})`
  i skończonej energii, `lim_{R→∞} ∮_{S_R} n_μ J^μ dS = 0`.
- Export:
  - struktura dowodu (kroki),
  - warunki techniczne (weighted Sobolev),
  - miejsca zależne od gauge-fixing.

## Wyjście
- `generated/p1661_s611_strict_h2_boundary_flux_lemma_scaffold_summary.json`
