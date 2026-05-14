# P1575 S525 Strict Global Error Propagation Bound Packet (No Legacy Bridge)

Status: `P1575_EXECUTED_GLOBAL_ERROR_PROPAGATION_BOUND_CANDIDATE`
As of: `2026-05-14`

## Cel

Po `P1574` formalizujemy `T1574B`:

- uporządkowanie geometrii overlap,
- estymacja globalnego boundu propagacji błędu przez cały łańcuch,
- wskazanie regionów dominujących błąd.

## Konstrukcja strict-only

1. Porządkujemy punkty overlap leksykograficznie `(omega,phi,beta,eta)`.
2. Dla każdego kroku wyznaczamy lokalny operator propagacji błędu
   `B_i = ||J_i^{-1}||_inf * g_i`, gdzie `g_i` to lokalny gain łańcucha
   `coefficients -> lagrangian -> eom`.
3. Globalny bound bierzemy jako `B_global = max_i B_i` i raportujemy top-k
   punktów dominujących.

## Kryterium PASS/FAIL

- `PASS_T1574B_CANDIDATE` jeśli `B_global` jest skończony i poniżej progu
  operacyjnego na całym overlap.

## Wynik

`FAIL_T1574B_CANDIDATE` (global bound still above operational threshold).

## Brakujące obiekty do final strict-core closure

1. `T1575A`: theorem odporności boundu na perturbacje EOM.
2. `W1575B`: witness zgodności boundu z celem ToE `F_nadsoliton => L_SM + L_GR`.
3. `T1575C`: końcowy theorem sklejenia globalnego.

## Rekomendacja następnego uczciwego kroku

Uruchomić `P1576`: perturbacyjny stress-test EOM z lokalnym dociążeniem punktów
dominujących i formalizacja `T1575A`.

## Omówienie dla laika

To licznik bezpieczeństwa błędu dla całego modelu: mówi, jak bardzo błąd na
wejściu może urosnąć po przejściu przez wszystkie etapy obliczeń.
