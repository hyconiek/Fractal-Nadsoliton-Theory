# P1571 S521 Strict Semiglobal Stability Scan Packet (No Legacy Bridge)

Status: `P1571_EXECUTED_SEMIGLOBAL_STABILITY_SCAN_T1570C_CANDIDATE`
As of: `2026-05-14`

## Cel

Wykonać kolejny uczciwy krok po `P1570`:

- semiglobalny skan stabilności inwersji na domenie admissible strict,
- wyznaczyć najgorszy punkt uwarunkowania,
- przygotować materiał pod `T1570C` (uniform stability bound).

## Konfiguracja strict-only

Mapa obserwowalna pozostaje:
`(lambda_sm_eff, kappa_gr_eff, epsilon_mix_eff, xi_phase_curv)`.

Domena robocza (siatka):
- `omega ∈ [0.16, 0.36]`,
- `phi ∈ [0.08, 0.24]`,
- `beta ∈ [0.85, 1.15]`,
- `eta ∈ [1.55, 1.95]`.

## Kryterium PASS/FAIL

- `PASS_T1570C_CANDIDATE`:
  1. determinant Jacobianu nie zmienia znaku na całej siatce,
  2. `min |det(J)| > det_threshold`,
  3. `max ||J^{-1}||_inf < cond_threshold`.

## Wynik

`FAIL_T1570C_CANDIDATE` (worst-point conditioning still above threshold).

## Brakujące obiekty do final strict-core closure

1. `T1571A`: theorem ciągłego patchingu kart lokalnych na pełnej domenie.
2. `W1571B`: witness zgodności boundów z pełnym `L_SM + L_GR` bundle.
3. `T1571C`: eksport globalnej granicy błędu dla klasy pomiarów EOM.

## Rekomendacja następnego uczciwego kroku

Uruchomić `P1572`: podział domeny (chart patching) z odcięciem obszaru złego
uwarunkowania i budową regionalnych boundów błędu inwersji.

## Omówienie dla laika

Sprawdziliśmy nie tylko jeden punkt, ale cały obszar parametrów. To jak test
samochodu nie tylko na prostym odcinku, ale na całej trasie: dzięki temu
wiemy, gdzie model jest najbardziej „wrażliwy”; tam trzeba zastosować podział
domeny zanim ogłosimy jednolity globalny bound.
