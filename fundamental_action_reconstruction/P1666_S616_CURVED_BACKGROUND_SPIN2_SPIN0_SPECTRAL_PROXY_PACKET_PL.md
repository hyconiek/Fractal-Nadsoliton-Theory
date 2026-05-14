# P1666 / S616 Curved-Background Spin-2/Spin-0 Spectral Proxy Packet

Status: `P1666_EXECUTED_CURVED_BACKGROUND_SPECTRAL_PROXY_NO_FALSE_PASS`
As of: `2026-05-14`

## Cel

Rozszerzyć `S615` o klasę teł zakrzywionych (FRW-proxy) i sprawdzić spektrum
spin-2/spin-0 na poziomie operacyjnym bez udawania theorem-level domknięcia.

## Tor strict-only

`K_strict -> coeff_GR/QG -> L_GR+QG_ct -> propagator poles/residues on {Minkowski, FRW_proxy}`

## Co testujemy

1. `m2_spin2(Rbg)` i `m2_spin0(Rbg)` dla dwóch teł,
2. brak tachionów: `m2 > 0`,
3. znak proxy-residuów spin-2/spin-0,
4. stabilność znaków między tłami.

## No-false-pass

Wynik eksportujemy jako `OPEN_OBLIGATION` nawet przy lokalnie dobrych liczbach,
bo nadal brak formalnego globalnego dowodu unitarności i renormalizacji.

## Rekomendowany następny uczciwy krok

`S617`: przejść z proxy do jawnej operatorowej diagonalizacji pełnego kernela kwadratowego i policzyć rzeczywiste residua (z projekcjami i-fix gauge) na klasie teł.

## Omówienie dla laika

To jak test auta na dwóch nawierzchniach: sucha i mokra.
Sprawdzamy, czy zachowanie trybów drgań nie staje się niefizyczne po zmianie tła,
ale pełna homologacja nadal wymaga głębszego dowodu matematycznego.
