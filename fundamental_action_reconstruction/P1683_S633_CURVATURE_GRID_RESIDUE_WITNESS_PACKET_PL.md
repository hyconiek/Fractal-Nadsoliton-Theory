# P1683 / S633 Curvature-Grid Residue Witness Packet (strict-only)

Status: `P1683_EXECUTED_WITNESS_EXPORT_NO_FALSE_PASS`
As of: `2026-05-14`

## Cel
Wykonać S633: wyeksportować witness dla `L_residue_curvature_uniform` poprzez
siatkę krzywizn i symboliczne ograniczenie znaku `det(K_spin2/0)` dla strict
`L_total`, bez bridge do legacy.

## Pełny tor fizyczny
`K_strict -> theta_strict -> coeff -> full L_total -> linearized EOM -> kinetic blocks -> residue witness`.

## Full Lagrangian anchor
`L_total = L_gauge + L_fermion + L_higgs + L_yukawa + L_gravity + L_mix`,
`L_gravity = (M_Pl^2/2)R + c1 R^2 + c2 R_{μν}R^{μν}`,
`L_mix = ξ H^† H R + strict counterterm set`.

## Zakres eksportu
- curvature grid: `R = {-R_max, -R_mid, 0, R_mid, R_max}`,
- testy znaków: `Z2(R)>0`, `Z0(R)>0`, `det(K_spin2)>0`, `det(K_spin0)>0` (proxy-level),
- status każdego punktu siatki + agregat globalny,
- jawne ograniczenie: to jeszcze witness eksportowy, nie theorem closure.

## Wynik
Eksport: `generated/p1683_s633_curvature_grid_residue_witness.json`.
Status globalny: `OPEN_OBLIGATION`.

## Braki do strict-core closure
1. przejście proxy-sign -> theorem-level inequality proof,
2. amplitude optical-theorem integration,
3. pełne sklejenie z background-independence cocycle theorem.

## Następny uczciwy krok
`S634`: podnieść witness do lemma-proof sketch z jawnie wyprowadzonymi
nierównościami funkcji `R` i warunkami na `(c1,c2,xi)`.

## Omówienie dla laika
To jak test odporności teorii na różne „krzywizny wszechświata”. Sprawdzamy,
czy podstawowe wielkości nie zmieniają znaku w niebezpieczny sposób. Wynik mówi,
że mamy mocniejszy dowód roboczy, ale jeszcze nie finalny dowód matematyczny.
