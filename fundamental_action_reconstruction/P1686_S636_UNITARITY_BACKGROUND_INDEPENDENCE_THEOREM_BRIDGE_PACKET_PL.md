# P1686 / S636 Unitarity + Background-Independence Theorem-Bridge Packet (strict-only)

Status: `P1686_EXECUTED_BRIDGE_SCAFFOLD_NO_FALSE_PASS`
As of: `2026-05-14`

## Cel
Wykonać S636: spiąć certyfikat granic S635 z globalnym łańcuchem cocycle
(background independence) i theorem-level statement unitarności, cały czas w
strict-only (`K_strict`) bez legacy bridge.

## Pełny tor fizyczny (przód + tył)
Przód:
`K_strict -> coeff(c1,c2,xi,...) -> full L_total -> EOM -> residue bounds -> unitarity statement -> cocycle/global BI bridge`.
Tył:
`global BI constraints -> overlap/cocycle consistency -> admissible unitarity domain -> coeff-sector -> strict-kernel consistency`.

## Full Lagrangian anchor
`L_total = L_gauge + L_fermion + L_higgs + L_yukawa + L_gravity + L_mix`,
`L_gravity = (M_Pl^2/2)R + c1 R^2 + c2 R_{μν}R^{μν}`,
`L_mix = ξ H^† H R + strict counterterm set`.

## Most theorem-level (scaffold)
- wejście A: `S635` bound certificate dla `(c1,c2,xi)` + boundary map,
- wejście B: `S629` cocycle scaffold (`L_overlap1/L_overlap2/L_cocycle3`),
- cel C: skomponować warunkowy theorem object:
  `T_unitarity_BI_bridge_strict : (A and B and lemma-set) => strict global consistency statement`.

## Wynik
Eksport: `generated/p1686_s636_unitarity_background_independence_theorem_bridge.json`.
Status: `OPEN_OBLIGATION`.

## Braki do final strict-core closure
1. brak pełnych dowodów `L_overlap1/L_overlap2/L_cocycle3`,
2. brak theorem-level optical-theorem amplitude proof,
3. brak globalnej pełnej inwersji `EOM <-> L_total` dla sektora sprzężonego SM+GR.

## Następny uczciwy krok
`S637`: uruchomić dowodowy pakiet dla `L_overlap1` i podłączyć go do
`T_unitarity_BI_bridge_strict` jako pierwszy realny zamykany blok.

## Omówienie dla laika
To etap łączenia dwóch dużych wymagań: teoria ma być jednocześnie unitarna
(i bez „duchów”) oraz niezależna od wyboru tła czasoprzestrzeni. Mamy już most
koncepcyjny i formalny szkielet, ale jeszcze nie komplet dowodów.
