# P1569 S519 Strict Local Uniqueness Theorem Candidate Packet (No Legacy Bridge)

Status: `P1569_EXECUTED_LOCAL_UNIQUENESS_THEOREM_CANDIDATE_EXPORT`
As of: `2026-05-14`

## Cel

Po `P1568` formalizujemy następny uczciwy krok:

`T1568D` — lokalna jednoznaczność inwersji mapy 4-obserwowalnej
`(lambda_sm_eff, kappa_gr_eff, epsilon_mix_eff, xi_phase_curv)`
na parametry strict `(omega, phi, beta, eta)`.

## Konstrukcja strict-only

Mapa obserwowalna:

- `lambda_sm_eff = |cos(phi)|/(1+0.3*beta)`
- `kappa_gr_eff = (omega^2 + 0.5*eta)/(0.2 + |cos(phi)|)`
- `epsilon_mix_eff = (beta*eta + 0.1*omega)/(1+|sin(phi)|)`
- `xi_phase_curv = omega*eta + sin(phi)`

Theorem-level candidate (lokalny, numeryczny):

1. `det(J4)` mapy 4->4 w punkcie strict jest oddalony od zera,
2. znak `det(J4)` jest stabilny na małej kuli roboczej,
3. minimalna wartość `|det(J4)|` na próbce lokalnej > próg bezpieczeństwa.

## Kryterium PASS/FAIL

- `PASS_T1568D_CANDIDATE`: wszystkie 3 warunki spełnione.
- `FAIL_T1568D_CANDIDATE`: co najmniej jeden warunek niespełniony.

## Wynik

`PASS_T1568D_CANDIDATE` — lokalna odwracalność mapy 4-obserwowalnej jest
numerycznie wsparta na zadanej domenie roboczej strict.

## Dalsze brakujące obiekty do final strict-core closure

1. `E1569A`: zamknięty analityczny eksport Jacobianu/Hessianu (bez różnic
   skończonych).
2. `T1569B`: formalny dowód ograniczenia błędu inwersji (stability envelope).
3. `T1569C`: globalny patching lokalnych kart inwersji do szerszej domeny ToE.

## Rekomendacja następnego uczciwego kroku

Uruchomić `P1570`: eksport `E1569A` (analityczny Jacobian/Hessian) i spinanie
z estimate'm stabilności `T1569B`.

## Omówienie dla laika

Dodany 4-ty sygnał działa jak brakujący wymiar informacji. Dzięki temu w
najbliższym otoczeniu punktu pracy modelu można już odróżniać jednoznacznie,
które parametry strict wygenerowały dany zestaw efektów fizycznych.
