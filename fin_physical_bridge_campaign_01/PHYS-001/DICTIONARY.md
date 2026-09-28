# PHYS-001 — exact Curie–Weiss dictionary

## Dokładna równoważność
Dla etykiet `sigma_a in Z_12` definiujemy wyłącznie jako słownik modelowy

`H(sigma)=-(J_c/(2N)) sum_{a,b} A7[sigma_a,sigma_b] - h sum_a 1[sigma_a=0]`.

Przy `beta_th=1/(k_B T)`, `g=beta_th J_c` i `vartheta=beta_th h=kappa/N` miara kanoniczna daje dla liczebności `n` dokładnie

`pi(n) ∝ N!/prod_j n_j! * exp[(g/(2N)) n^T A7 n + vartheta n_0]`.

Współczynnik `1/2` usuwa podwójne liczenie par w sumie po `a,b`; `1/N` jest skalowaniem mean-field zapewniającym energię ekstensywną. Czynnik multinomialny jest dokładną liczbą mikrostanów o danych liczebnościach.

Po usunięciu aktualizowanej kopii `i`, `m=n-e_i`. Ponieważ A7 ma stałą przekątną, składnik samodiagonalny nie zależy od proponowanej etykiety i znika w normalizacji:

`q_j(m)=softmax_j[(g/N)(A7 m)_j + vartheta 1[j=0]]`.

Stąd count-rate dla heat-bath wynosi `nu n_i q_j(m)`; `nu=1` w repo jest wyborem jednostki generator-time, a nie sekundą.

## Generalized Curie–Weiss, Potts i vector-spin
To jest dokładny macierzowy/generalized q-state Curie–Weiss. Nie jest automatycznie najprostszym Potts Hamiltonianem `-J delta_{sigma_a,sigma_b}`, bo A7 nie ma postaci jednego centered Potts projector. Ponieważ `A7=X7 X7^T`, można pisać `xi_j=X7[j,:]` i

`H=-(J_c/(2N)) |sum_a xi_{sigma_a}|^2 - h n_0`.

Jest to więc także dyskretny wektorowy model mean-field z 12 dozwolonymi wektorami wewnętrznymi w `R^7`. `R^7` jest przestrzenią cech/spinów tego modelu, nie wyprowadzoną 7-wymiarową przestrzenią fizyczną.

## Gauge skali
Transformacja `A -> c A`, `J_c -> J_c/c` pozostawia Hamiltonian parowy bez zmian. Zatem sama skala A nie jest identyfikowalna bez konwencji dla `J_c`/`g`; kształt widma po normalizacji jest osobną informacją.

## Granica roszczenia
Równoważność miary Gibbsa jest matematyczna (M). Wybór heat-bath jest dodatkową regułą kinetyczną. Utożsamienie `J_c,h,T,nu` z wielkościami konkretnego urządzenia wymaga osobnego mostu fizycznego/kalibracji; target Hamiltonian nie jest przez to energią elektryczną sprzętu.
