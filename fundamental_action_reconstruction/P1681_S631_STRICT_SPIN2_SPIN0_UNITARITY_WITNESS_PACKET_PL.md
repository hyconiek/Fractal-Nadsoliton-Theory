# P1681 / S631 Strict Spin-2/Spin-0 Unitarity Witness Packet (strict-only)

Status: `P1681_EXECUTED_WITNESS_SCAFFOLD_NO_FALSE_PASS`
As of: `2026-05-14`

## Cel
Wykonać krok S631 wskazany w P1680: zbudować theorem-level scaffold dla
unitarity witness sektora spin-2/spin-0 na curved background, ściśle w torze
`K_strict -> coeff -> L_total -> EOM`, bez bridge do legacy.

## Pełny tor fizyczny (jawny)
1. `K_strict(d)` ustala strict-parametry `theta_strict`.
2. `theta_strict` generuje współczynniki `{Z2,Z0,m2_sq,m0_sq,c1,c2,xi,...}`.
3. Współczynniki osadzane są w pełnym `L_total = L_SM + L_GR + L_mix`.
4. Z `L_total` wyprowadzamy EOM i liniaryzujemy wokół curved background.
5. Z operatora kwadratowego konstruujemy projekcję spin-2/spin-0 i testy unitarności.

## Full Lagrangian anchor (nieszkieletowy)
`L_total = L_gauge + L_fermion + L_higgs + L_yukawa + L_gravity + L_mix`,
z `L_gravity = (M_Pl^2/2)R + c1 R^2 + c2 R_{μν}R^{μν}` oraz
`L_mix = ξ H^† H R + strict counterterm set`.

## Eksport świadka (scaffold)
Eksport: `generated/p1681_s631_strict_spin2_spin0_unitarity_witness_scaffold.json`.
Zawiera:
- mapę parametrów strict -> sektor unitarności,
- warunki konieczne (`Z2>0`, `Z0>0`, `m2_sq>=0`, `m0_sq>=0`),
- status theorem-level dowodu (nadal otwarty),
- jawne luki (residuum globalne + optical theorem + background-globality).

## Klasyfikacja uczciwości
`OPEN_OBLIGATION` — nie deklarujemy domknięcia unitarności ani ToE.

## Braki do final strict-core closure
1. global proof dodatniości residuów na pełnej klasie backgroundów,
2. optical-theorem compatibility dla amplitud strict EFT z miksami SM-GR,
3. zapięcie z renormalization/background-independence theorem chain,
4. global Helmholtz `EOM -> L_total` theorem spin-coupled.

## Następny uczciwy krok
`S632`: wyeksportować lemma-level package dla (1): curvature-uniform residue bound
w atlasie niecyklicznym kompatybilnym z `F3 kernel-split-robust`.

## Omówienie dla laika
Sprawdzamy, czy „cząstki” wynikające z teorii mają sens fizyczny: czy nie pojawiają się
niestabilne duchy (ghosty) i czy rachunek prawdopodobieństwa oddziaływań jest spójny.
Na razie mamy rygorystyczny plan i warunki konieczne, ale jeszcze nie pełny dowód.
