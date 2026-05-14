# P1660 / S610 — Strict H2 gauge-scalar boundary-condition candidate export

## Cel
Po P1659: wyeksportować jawny kandydat warunku brzegowego, który zeruje
`∫_∂Ω n_μ J^μ` dla klasy pól o skończonej energii w bloku gauge-scalar.

## Zakres
- strict-only, bez legacy bridge,
- wynik `PARTIAL` (kandydat + obowiązki dowodowe),
- brak deklaracji domknięcia H2/H3/H4 i bram QG.

## Konstrukcja
- Kandydat klasy pól:
  - `A_μ = O(r^{-1-ε})`,
  - `D_μ φ = O(r^{-2-ε})`,
  - lokalna skończoność energii i ładunku.
- Wniosek roboczy: wkład brzegowy maleje jak `r^{-1-ε}` i zanika dla `r -> ∞`.

## Wyjście
- `generated/p1660_s610_strict_h2_gauge_scalar_boundary_condition_candidate_summary.json`
