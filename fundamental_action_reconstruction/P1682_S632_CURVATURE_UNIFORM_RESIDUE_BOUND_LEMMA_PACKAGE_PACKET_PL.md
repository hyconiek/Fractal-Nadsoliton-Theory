# P1682 / S632 Curvature-Uniform Residue-Bound Lemma Package (strict-only)

Status: `P1682_EXECUTED_LEMMA_PACKAGE_SCAFFOLD_NO_FALSE_PASS`
As of: `2026-05-14`

## Cel
Zrealizować krok S632: przygotować lemma-level pakiet dla globalnej dodatniości
residuów spin-2/spin-0 na niecyklicznym atlasie backgroundów (klasa
`kernel-split-robust`), w torze strict-only bez legacy bridge.

## Łańcuch fizyczny (dwukierunkowy anchor)
`K_strict -> coeff(theta) -> full L_total -> linearized EOM -> residue lemmas -> unitarity theorem chain`.
Wsteczny kierunek audytowy: `EOM constraints -> admissible coeff-sector -> consistency with K_strict parameter map`.

## Full Lagrangian anchor
`L_total = L_gauge + L_fermion + L_higgs + L_yukawa + L_gravity + L_mix`,
`L_gravity = (M_Pl^2/2)R + c1 R^2 + c2 R_{μν}R^{μν}`,
`L_mix = ξ H^† H R + strict counterterm set`.

## Pakiet lematów (S632 scaffold)
- `L_residue_local_pos`: lokalna dodatniość residuów dla chartu `U_i`.
- `L_residue_overlap_stability`: zgodność znaków residuów na overlapach `U_i∩U_j`.
- `L_residue_curvature_uniform`: jednorodna granica w klasie krzywizn `|R|<=R_max`.
- `L_residue_global_glue`: przejście lokalne -> globalne w atlasie niecyklicznym.

## Wynik
Eksport: `generated/p1682_s632_curvature_uniform_residue_bound_lemma_package.json`.
Status: `OPEN_OBLIGATION` (brak twierdzenia finalnego).

## Braki do strict-core closure
1. formalny dowód `L_residue_curvature_uniform` bez luki granicznej,
2. sprzężenie z optical-theorem amplitude constraints,
3. pełne sklejenie z renormalization + background-independence theorem chain.

## Następny uczciwy krok
`S633`: wygenerować dowodowy witness dla `L_residue_curvature_uniform` na siatce
krzywizn + symbolicznym ograniczeniu sign(det kinetic block) dla spin-2/spin-0.

## Omówienie dla laika
To etap „sprawdzania stabilności” teorii w różnych warunkach grawitacyjnych.
Chcemy mieć pewność, że teoria pozostaje fizyczna nie tylko lokalnie, ale globalnie,
nawet gdy tło czasoprzestrzeni się zmienia.
