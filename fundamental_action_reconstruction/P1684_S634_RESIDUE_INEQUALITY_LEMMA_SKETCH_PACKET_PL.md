# P1684 / S634 Residue-Inequality Lemma Sketch Packet (strict-only)

Status: `P1684_EXECUTED_LEMMA_SKETCH_NO_FALSE_PASS`
As of: `2026-05-14`

## Cel
Wykonać S634: podnieść witness S633 do szkicu lemma-proof z jawnymi
nierównościami w krzywiźnie `R` i warunkami na `(c1,c2,xi)` dla sektora
spin-2/spin-0, bez bridge do legacy.

## Pełny tor fizyczny
`K_strict -> theta_strict -> coeff(c1,c2,xi,...) -> full L_total -> linearized EOM -> residue inequalities`.
Kierunek wsteczny audytowy: `ineq-residue constraints -> admissible coeff domain -> consistency with strict kernel map`.

## Full Lagrangian anchor
`L_total = L_gauge + L_fermion + L_higgs + L_yukawa + L_gravity + L_mix`,
`L_gravity = (M_Pl^2/2)R + c1 R^2 + c2 R_{μν}R^{μν}`,
`L_mix = ξ H^† H R + strict counterterm set`.

## Szkic nierówności (lemma-level)
- `Z2(R) = Z2_0 + a2*R >= z2_min > 0`,
- `Z0(R) = Z0_0 + a0*R >= z0_min > 0`,
- `det(K2)(R) >= k2_min > 0`,
- `det(K0)(R) >= k0_min > 0`,
na domenie `|R|<=R_max`.

## Wynik
Eksport: `generated/p1684_s634_residue_inequality_lemma_sketch.json`.
Status: `OPEN_OBLIGATION` (proof sketch, nie theorem closure).

## Braki do strict-core closure
1. formalizacja dowodu nierówności z pełnym rachunkiem wariacyjnym,
2. spięcie z optical theorem dla amplitud,
3. pełna integracja z renormalization + background independence.

## Następny uczciwy krok
`S635`: wygenerować symboliczny pakiet granic (`bound certificate`) dla
`(c1,c2,xi)` wraz z kontrprzykładami granicznymi i mapą domeny dopuszczalnej.

## Omówienie dla laika
To przejście od „testu punktowego” do „reguł matematycznych” mówiących, kiedy
teoria zachowuje stabilność. Nadal nie zamykamy teorii na siłę, ale precyzujemy
warunki, które muszą być spełnione.
