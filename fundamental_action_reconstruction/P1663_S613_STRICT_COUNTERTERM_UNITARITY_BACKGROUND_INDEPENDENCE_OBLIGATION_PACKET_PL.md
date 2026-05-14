# P1663 / S613 Strict QG Obligation Packet: Counterterms + Unitarity + Background Independence

Status: `P1663_EXECUTED_STRICT_QG_OBLIGATION_PACKET_NO_FALSE_PASS`
As of: `2026-05-14`

## Cel

Następny uczciwy krok strict-only po `P1662`: przejść z deklaracji do jawnego, wykonywalnego eksportu
trzech obowiązków QG dla toru `K_strict -> coefficients -> L_full -> EOM`.

## Co dodajemy

1. **Counterterm witness (1-loop proxy)**
   - jawna baza operatorów: `R^2`, `R_{mu nu}R^{mu nu}`, `R_{mu nu rho sigma}R^{mu nu rho sigma}`,
   - jawny wektor współczynników kontrterminów `(delta_cR2, delta_cRic2, delta_cRiem2)`,
   - reguła strict-only mapowania z parametrów kernela.

2. **Unitarity witness (spin-2 proxy)**
   - warunek bez-duchowy (operacyjny): `Mpl2 > 0` oraz `cRic2 >= 0`,
   - raport `PASS/OPEN` bez fałszywego domknięcia theorem-level.

3. **Background-independence witness (variational proxy)**
   - zapis równań Eulera-Lagrange'a z tego samego `L_total` na dwóch tłach roboczych,
   - wymóg: zgodna forma operatorowa po podstawieniu tła (`Minkowski`, `FRW_proxy`).

## Pełny łańcuch fizyczny (strict-only)

- `K_strict(omega,phi,beta,eta,A)`
  -> mapa współczynników EFT/GR,
- współczynniki
  -> pełny ansatz `L_SM + L_GR + L_mix + L_QG_ct`,
- `L_total`
  -> EOM sektorów `h`, `A_mu`, oraz macierz zobowiązań QG,
- EOM + constraints
  -> odwrotna kontrola identyfikowalności współczynników.

## Zasada no-false-pass

`P1663` może dać:

- `PASS` tylko dla lokalnych testów algebraiczno-numerycznych,
- `OPEN_OBLIGATION` dla nierozstrzygniętych theorem-level aspektów renormalizacji/unarności/background independence.

## Wynik

Skrypt `p1663_strict_qg_obligation_witness.py` eksportuje raport JSON z:

- pełnym zestawem współczynników strict,
- bazą kontrterminów i liczbami proxy,
- testem unitarności proxy,
- dwoma evaluacjami tła i residualami EOM.

## Rekomendowany następny uczciwy krok

`S614`: zastąpić proxy-kontrterminy jawnie policzonymi strukturami z rozwinięcia 1-loop (heat-kernel / background-field) i dodać operatorowy test unitarności propagatora spin-2.

## Omówienie dla laika

To jak sprawdzanie nowego samolotu:

- najpierw budujesz kompletny model (lagranżian),
- potem liczysz, jak powinien latać (równania ruchu),
- na końcu robisz testy bezpieczeństwa (czy nie ma niestabilności i czy działa w różnych warunkach).

W tym kroku te testy są już uruchomione i mierzalne, ale jeszcze nie są końcową certyfikacją fizyki kwantowej grawitacji.
