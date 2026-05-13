# P1402 — STRICT-ONLY PC4 PROVIDER BASELINE DESIGN (PL)

## Cel kroku
Po aktywacji `PC4` w `P1401` zamrażamy nową bazę provider-class dla toru strict-only
`F_Nadsoliton ⇒ L_SM + L_GR`, bez bridge do legacy i bez powrotu do pętli `PC3`.

## Zakres strict-only
- `legacy_bridge_used = false`
- `inherits_from_pc3_loop = false`
- `noncyclic_provider_shift = enforced`

## Baseline PC4 (v1)
- `provider_class_id = PC4_strict_phase_locked_selector_anchor_v1`
- `anchor_family = A4_noncyclic_phase_locked_stabilizers`
- `epsilon_sign_v1 = 0.05` (zamrożone)
- `epsilon_drift_v1 = 0.04` (zamrożone)
- `first_run_contract = dual_metric_required`

## Status kroku
- `PC4_BASELINE_STATUS := FROZEN_READY_FOR_FIRST_RUN`
- `L_B1_03_EXPORT_STATUS := NOT_EXPORTED`
- `B1_STATUS := OPEN`

## Decyzja profesorska (następny uczciwy krok)
Uruchomić `P1403_STRICT_ONLY_PC4_FIRST_DUAL_METRIC_RUN` na zamrożonych progach i wydać jawny PASS/FAIL bez retuningu kryteriów po fakcie.

## Omówienie dla laika
To start nowej wersji modelu: ustawiamy zasady i progi zanim zaczniemy test. Dzięki temu wynik będzie uczciwy i porównywalny, a nie „dopasowany po fakcie”.
