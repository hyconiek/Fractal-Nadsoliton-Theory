# P1664 / S614 Strict Full Lagrangian Manifest + Bidirectional Inversion Packet

Status: `P1664_EXECUTED_FULL_LAGRANGIAN_MANIFEST_NO_FALSE_PASS`
As of: `2026-05-14`

## Cel

Dostarczyć **pełny manifest lagranżianu** (nie-szkieletowy) w strict-core oraz wykonać
obustronny test toru:

`K_strict -> coefficients -> L_full -> EOM-constraints -> recovered kernel`.

## Pełny lagranżian strict-core (manifest)

`L_total = L_GR + L_gauge + L_H + L_fermion + L_Yukawa + L_mix + L_QG_ct`

z sektorami:

- `L_GR = sqrt(-g)[ (Mpl2/2)R + cR2 R^2 + cRic2 R_{mu nu}R^{mu nu} + cRiem2 R_{mu nu rho sigma}R^{mu nu rho sigma} ]`
- `L_gauge = -sqrt(-g)/4 * [ Z3 G^a_{mu nu}G^{a mu nu} + Z2 W^i_{mu nu}W^{i mu nu} + Z1 B_{mu nu}B^{mu nu} ]`
- `L_H = sqrt(-g)[ (D_mu H)^†(D^mu H) - muH2 H^†H - lambdaH(H^†H)^2 ]`
- `L_fermion = sqrt(-g) * sum_f i psi_f_bar gamma^a e_a^mu D_mu psi_f`
- `L_Yukawa = -sqrt(-g)[ y_u Q_bar H_tilde u_R + y_d Q_bar H d_R + y_e L_bar H e_R + h.c. ]`
- `L_mix = sqrt(-g)[ xiHR (H^†H)R + chiRG R*(G^2+W^2+B^2) ]`
- `L_QG_ct = sqrt(-g)[ delta_cR2 R^2 + delta_cRic2 Ricci^2 + delta_cRiem2 Riemann^2 ]`

## Obustronny test

Skrypt `p1664_strict_full_lagrangian_manifest_and_inversion.py`:

1. liczy mapę strict kernel -> coefficients,
2. buduje pełny manifest `L_total`,
3. wykonuje odwrócenie `coefficients -> kernel` dla `(beta, A, omega, eta)`,
4. raportuje błędy rekonstrukcji.

PASS lokalny tylko gdy błędy < tolerancji numerycznej.

## Wynik i rygor

- obustronny test przechodzi lokalnie (numerycznie),
- brak twierdzenia globalnego => status całej ścieżki nadal `OPEN_OBLIGATION` dla QG.

## Rekomendowany następny uczciwy krok

`S615`: wyprowadzić jawny operatorowy propagator spin-2 ze składników `R^2/Ricci^2/Riemann^2` i dopisać test resztowy unitarności (brak biegunów o ujemnej normie) na reprezentatywnej klasie teł.

## Omówienie dla laika

Tutaj mamy pełną „listę części i połączeń” teorii (pełny lagranżian), a nie tylko szkic.
Następnie sprawdzamy, czy z obserwowanego działania modelu da się odzyskać ustawienia startowe.
To ważne, bo pokazuje że model działa w obie strony: od założeń do przewidywań i z powrotem.
