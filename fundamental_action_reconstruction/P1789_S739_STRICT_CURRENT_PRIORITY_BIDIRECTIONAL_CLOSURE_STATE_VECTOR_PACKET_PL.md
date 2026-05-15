# P1789 S739 Strict Current Priority Bidirectional Closure State Vector Packet (PL)

## Cel

Skonsolidować jeden stan wektorowy dla całego bieżącego priorytetu
(Forward + Reverse + Theorem Gates), bez nadmiarowego rozmnażania checkpointów.

## Scope

Forward lane:
- `K_strict -> coefficients -> full non-skeleton L_total -> covariant nonproxy EOM bundle`.

Reverse lane:
- `H1` weak-form consistency,
- `EL_g - E_{μν}` componentwise residual,
- promotion policy do theorem-level gates.

## Current state vector (strict-only, no false-pass)

1. `SV1_E_A_mu_nonproxy_covariant_explicit` -> `OPEN_COMPONENTWISE_REQUIRED`.
2. `SV2_E_H_nonproxy_covariant_explicit` -> `OPEN_COMPONENTWISE_REQUIRED`.
3. `SV3_EL_g_nonproxy_explicit` -> `OPEN_COMPONENTWISE_REQUIRED`.
4. `SV4_boundary_term_control` -> `PARTIAL_OPEN`.
5. `SV5_H1_4D_weak_form` -> `OPEN`.
6. `SV6_Bianchi_Ward_global` -> `OPEN_LOCKED`.
7. `SV7_BRST_global` -> `OPEN_LOCKED`.
8. `SV8_Cutkosky_global` -> `OPEN_LOCKED`.

## Discipline of interpretation

- `LOCAL PASS` nigdy nie implikuje `GLOBAL PASS`.
- `SCAFFOLD/TEMPLATE` nigdy nie implikuje `FULL EXPORT`.
- Każdy brak jawnego residual/witness pozostaje `OPEN`.

## Minimal bidirectional closure condition

Warunek minimalny na uczciwe przejście do rozmowy o theorem-level promotion:

- jednoczesny komponentowy run `E_A^μ`, `E_H`, `EL_g` na wspólnej rodzinie teł,
- jawny ledger `H1` oraz `EL_g-E_{μν}`,
- dopiero potem sekwencja `BW -> BRST -> Cutkosky`.

## Ryzyka false-pass

1. Promocja globalna po pojedynczym lokalnym sukcesie.
2. Pominięcie boundary-control przy interpretacji H1.
3. Przeskok do BRST/Cutkosky bez domknięcia BW.

## Następny uczciwy krok

Wykonać jeden skoordynowany run-pack dla `SV1..SV5` i wydać wynik wyłącznie
jako `PASS_ZERO` albo `OPEN_OBSTRUCTION_WITH_TRACE`; bez statusów pośrednio-sugerujących domknięcie globalne.

## Objaśnienie dla laika

To tablica kontrolna projektu: pokazuje dokładnie, które lampki są jeszcze czerwone,
a które dopiero żółte — bez udawania, że cel już osiągnięto.
