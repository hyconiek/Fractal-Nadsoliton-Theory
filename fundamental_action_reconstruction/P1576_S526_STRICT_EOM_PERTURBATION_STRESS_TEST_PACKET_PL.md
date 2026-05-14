# P1576 S526 Strict EOM Perturbation Stress-Test Packet (No Legacy Bridge)

Status: `P1576_EXECUTED_EOM_PERTURBATION_ROBUSTNESS_CANDIDATE_T1575A`
As of: `2026-05-14`

## Cel

Po `P1575` sprawdzamy odporność globalnego boundu na perturbacje EOM
w punktach dominujących błąd.

## Konstrukcja strict-only

1. Bierzemy top punkty dominujące z `P1575`.
2. Nakładamy perturbacje EOM o amplitudach `1e-4, 3e-4, 1e-3`.
3. Mierzymy względny wzrost lokalnego `B_local`.

## Kryterium PASS/FAIL

- `PASS_T1575A_CANDIDATE` jeśli maksymalny względny wzrost `B_local`
  pozostaje poniżej progu odporności.

## Wynik

`FAIL_T1575A_CANDIDATE` (relative growth above robustness threshold).

## Brakujące obiekty do final strict-core closure

1. `W1576A`: witness zgodności perturbacyjnej odporności z full `L_SM + L_GR`.
2. `T1576B`: theorem globalnego sklejenia z odpornością szumową.
3. `T1576C`: final strict-core closure theorem bundle.

## Rekomendacja następnego uczciwego kroku

Uruchomić `P1577`: eksport `W1576A` i integracja perturbacyjnej odporności z
pełnym bundlem ToE.

## Omówienie dla laika

To test „co jeśli pomiar jest zaszumiony”: sprawdzamy, czy drobne błędy w
równaniach ruchu nie rozwalają całej rekonstrukcji parametrów.
