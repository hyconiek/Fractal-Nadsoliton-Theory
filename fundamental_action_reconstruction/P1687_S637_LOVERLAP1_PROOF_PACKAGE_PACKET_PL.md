# P1687 / S637 L_overlap1 Proof-Package Packet (strict-only)

Status: `P1687_EXECUTED_PROOF_PACKAGE_SCAFFOLD_NO_FALSE_PASS`
As of: `2026-05-14`

## Cel
Wykonać S637: uruchomić pakiet dowodowy dla `L_overlap1` i podłączyć wynik do
`T_unitarity_BI_bridge_strict`, zachowując strict-only tor bez legacy bridge.

## Pełny tor fizyczny
`K_strict -> coeff -> full L_total -> EOM -> local overlap operators -> L_overlap1 proof package -> bridge attachment`.

## Full Lagrangian anchor
`L_total = L_gauge + L_fermion + L_higgs + L_yukawa + L_gravity + L_mix`,
`L_gravity = (M_Pl^2/2)R + c1 R^2 + c2 R_{μν}R^{μν}`,
`L_mix = ξ H^† H R + strict counterterm set`.

## Zakres pakietu dowodowego L_overlap1
- definicja operatorów przejścia `T_01`, `T_12` na overlapach `U0∩U1`, `U1∩U2`,
- warunek zgodności: `||T_12∘T_01 - T_02|| <= eps_overlap`,
- jawna lista założeń lokalnej regularności i domeny funkcji,
- status podłączenia do obiektu `T_unitarity_BI_bridge_strict`.

## Wynik
Eksport: `generated/p1687_s637_loverlap1_proof_package.json`.
Status globalny: `OPEN_OBLIGATION`.

## Braki do strict-core closure
1. formalny dowód nierówności overlap bez parametrycznego luzu `eps_overlap`,
2. propagacja wyniku do `L_overlap2` i `L_cocycle3`,
3. finalne sklejenie z optical theorem + renormalization constraints.

## Następny uczciwy krok
`S638`: wygenerować proof package dla `L_overlap2` z kompatybilnością
background-shift covariance na tych samych overlapach.

## Omówienie dla laika
To krok, w którym sprawdzamy, czy lokalne „przepisy” teorii są zgodne na
obszarach styku. Jeśli nie działają na styku, cała globalna teoria nie może być
spójna — dlatego ten etap jest krytyczny.
