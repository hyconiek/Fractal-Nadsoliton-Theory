# P1568 S518 Strict Underdetermined Orbit Witness And 4th Observable Packet (No Legacy Bridge)

Status: `P1568_EXECUTED_UNDERDETERMINED_ORBIT_WITNESS_AND_4TH_OBSERVABLE_CANDIDATE`
As of: `2026-05-14`

## Cel

Domknąć następny uczciwy krok po `P1567` na torze strict-only:

`kernel strict -> współczynniki -> Lagrangian -> EOM -> inwersja`.

Budujemy dwa brakujące obiekty:

1. `W1567A`: witness orbity niedookreślenia (różne parametry strict dają ten
   sam bundle obserwowalnych EOM).
2. `W1567B`: kandydat minimalnego 4-go obserwowalnego, który rozcina tę orbitę.

## Konstrukcja fizyczna (strict-only)

Nie używamy żadnego legacy bridge. Pracujemy lokalnie wokół punktu
`(omega, phi, beta, eta) = (0.18575, 0.16250, 1.0, 1.8)`.

- Bundle EOM-obserwowalnych: `(lambda_sm_eff, kappa_gr_eff, epsilon_mix_eff)`.
- Kandydat 4-go obserwowalnego: `xi_phase_curv = omega*eta + sin(phi)`
  (fazowo-krzywiznowy sygnał strict; niezależny od samego bundle 3x1).

## Kryterium PASS/FAIL

- `PASS_W1567A`: znaleziono co najmniej 2 różne zestawy parametrów strict,
  które odtwarzają ten sam bundle 3-obserwowalny w tolerancji numerycznej.
- `PASS_W1567B`: te same zestawy dają rozróżnialne wartości `xi_phase_curv`.

## Wynik

`PASS_W1567A_AND_PASS_W1567B_CANDIDATE`.

Wniosek: 3-obserwowalny bundle nie domyka pełnej inwersji; 4-ty obserwowalny
jest realnym kandydatem do domknięcia identyfikowalności strict-core.

## Brakujące dalsze obiekty do final strict-core closure

1. `E1568C`: analityczny eksport gradientu i Hessianu dla `xi_phase_curv`.
2. `T1568D`: theorem lokalnej jednoznaczności inwersji dla 4-obserwowalnego
   mapowania na admissible strict domain.
3. `T1568E`: theorem stabilności (warunek Lipschitz/condition number envelope).

## Rekomendacja następnego uczciwego kroku

Uruchomić `P1569`: formalizacja `T1568D` (lokalna jednoznaczność) z jawnym
zakresem domeny i warunkami brzegowymi dla parametrów strict.

## Omówienie dla laika

Pokazaliśmy, że sam zestaw 3 sygnałów z równań ruchu to za mało — różne
„wnętrza” modelu mogą wyglądać tak samo na zewnątrz. Dodanie czwartego,
niezależnego sygnału działa jak dodatkowy „odcisk palca”, który rozróżnia te
ukryte przypadki.
