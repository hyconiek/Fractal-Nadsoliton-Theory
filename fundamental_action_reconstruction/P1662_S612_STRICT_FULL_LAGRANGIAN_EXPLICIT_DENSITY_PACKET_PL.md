# P1662 / S612 — Strict full Lagrangian explicit density export (non-skeleton)

## Cel
Wyeksportować jedną jawnie rozwiniętą postać `L_total` dla toru
`K_strict -> współczynniki -> L_total -> EOM`, aby utrzymać pełny (nieszkieletowy) zapis
w aktualnym ciągu checkpointów strict-only.

## Zakres
- strict-only (`legacy_bridge_used=false`),
- bez claimu o final closure,
- obiekt ma służyć jako referencja do dalszych witnessów EOM<->L.

## Konstrukcja
`L_total = L_scalar + L_gauge + L_fermion + L_higgs + L_gravity + L_mix`, gdzie każdy składnik
jest zapisany jawnie z operatorami i sprzężeniami.

## Wyjście
- `generated/p1662_s612_strict_full_lagrangian_explicit_density_summary.json`
