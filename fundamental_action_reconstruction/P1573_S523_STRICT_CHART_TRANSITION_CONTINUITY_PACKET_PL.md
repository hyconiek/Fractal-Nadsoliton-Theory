# P1573 S523 Strict Chart-Transition Continuity Packet (No Legacy Bridge)

Status: `P1573_EXECUTED_CHART_TRANSITION_CONTINUITY_CANDIDATE`
As of: `2026-05-14`

## Cel

Po `P1572` formalizujemy `T1572B`:

- jawna mapa przejścia `good -> buffer` oraz `buffer -> bad`,
- test ciągłości na punktach overlap,
- eksport błędu sklejenia kart (gluing error).

## Konstrukcja strict-only

Dla każdego punktu `p=(omega,phi,beta,eta)` mamy lokalny deskryptor
`x = (lambda_sm_eff, kappa_gr_eff, epsilon_mix_eff, xi_phase_curv)`.

Mapy przejścia modelujemy lokalnie afinicznie:
`x_B = A_GB x_G + b_GB`, `x_D = A_BD x_B + b_BD`
(z kalibracją na wspólnych punktach overlap).

## Kryterium PASS/FAIL

- `PASS_T1572B_CANDIDATE` jeśli:
  1. overlap `good-buffer` i `buffer-bad` są niepuste,
  2. błąd sklejenia `L2` na overlapach < próg,
  3. kompozycja `good->buffer->bad` nie łamie ciągłości w granicy.

## Wynik

`PASS_T1572B_CANDIDATE`.

## Brakujące obiekty do final strict-core closure

1. `W1573A`: witness kompatybilności continuity z pełnym `L_SM + L_GR` bundle.
2. `T1573B`: theorem globalnego błędu po sklejeniu kart.
3. `T1573C`: theorem odporności patchingu na perturbacje EOM.

## Rekomendacja następnego uczciwego kroku

Uruchomić `P1574`: integracja continuity witness z pełnym łańcuchem
`kernel -> coefficients -> lagrangian -> EOM` i eksport `W1573A`.

## Omówienie dla laika

To etap „zszywania map”: upewniamy się, że gdy przechodzimy między strefami,
model nie robi skoku ani pęknięcia. Dzięki temu całość działa jak jedna,
spójna nawigacja, a nie zestaw luźnych kawałków.
