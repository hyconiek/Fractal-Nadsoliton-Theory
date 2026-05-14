# P1671 / S621 Dense Grid Variational Sensitivity + Bidirectional Proxy Packet

Status: `P1671_EXECUTED_DENSE_GRID_SENSITIVITY_AND_REVERSE_PROXY_NO_FALSE_PASS`
As of: `2026-05-14`

## Cel

Wykonać następny uczciwy krok po `P1670/S620` na torze strict-only:

`K_strict -> współczynniki -> pełny L_total -> równania ruchu`

z gęstszą siatką parametrów i jawnie policzoną wrażliwością residuów EOM na
`cR2,cRic2,cRiem2` oraz proxy-odwrotnością `EOM -> coeff`.

## Zakres rygoru

- strict-only (`legacy_bridge_used=false`),
- bez claimu o final closure,
- status pozostaje `OPEN_OBLIGATION` dopóki nie ma theorem-level dowodów dla:
  renormalizacji, unitarności spin-2/spin-0 i background independence.

## Co eksportujemy

1. Gęsta siatka 30 punktów `(beta,eta,omega,A,R)`.
2. Jawny pełny lagranżian sektorowy (nieszkieletowy) jako referencja operatorowa:
   `L_total = L_scalar + L_gauge + L_fermion + L_higgs + L_gravity + L_mix`.
3. Residua `metric/scalar/gauge` i lokalne pochodne numeryczne
   `d(residual)/d(cR2,cRic2,cRiem2)`.
4. Reverse-proxy: lokalna rekonstrukcja `alphaQ` z równania metrycznego
   (`EOM -> coeff`) oraz błąd rekonstrukcji.

## Wynik

Eksport: `generated/p1671_s621_dense_grid_variational_sensitivity_and_bidirectional_proxy.json`.

Interpretacja: forward chain jest obliczeniowo jawny i reprodukowalny, reverse chain
ma lokalną rekonstrukcję proxy, ale globalna odwracalność/theorem-level bridge
pozostają otwarte.

## Następny uczciwy krok

`S622`: zastąpić lokalny reverse-proxy theorem-level witness-em dla
`EOM -> L_total` (Helmholtz H2/H3/H4 na pełnej bazie operatorów strict) i spiąć
z dowodami QG (`renormalization/unitarity/background independence`).

## Omówienie dla laika

Zrobiliśmy gęstszy „test wytrzymałości” pełnego modelu strict:
nie tylko policzyliśmy równania ruchu, ale też sprawdziliśmy jak mocno wynik
zmienia się przy małej zmianie kluczowych współczynników grawitacyjnych.
Dodatkowo wykonaliśmy pierwszy krok „wstecz” (z równań do współczynników),
ale to nadal wersja robocza — pełny dowód matematyczny jeszcze przed nami.
