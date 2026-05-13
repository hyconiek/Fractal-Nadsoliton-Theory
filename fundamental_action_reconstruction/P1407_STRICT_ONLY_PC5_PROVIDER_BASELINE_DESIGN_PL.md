# P1407 — STRICT-ONLY PC5 PROVIDER BASELINE DESIGN (PL)

## Cel kroku
Po aktywacji `PC5` w `P1406` zamrażamy nową bazę provider-class dla strict-only toru
`F_Nadsoliton ⇒ L_SM + L_GR`, bez bridge do legacy.

## Zakres strict-only
- `legacy_bridge_used = false`
- `inherits_from_pc4_loop = false`
- `noncyclic_provider_shift = enforced`

## Baseline PC5 (v1)
- `provider_class_id = PC5_strict_dual_anchor_selector_v1`
- `anchor_family = A5_noncyclic_dual_anchor_stabilizers`
- `epsilon_sign_v1 = 0.05`
- `epsilon_drift_v1 = 0.04`
- `first_run_contract = dual_metric_required`

## Status kroku
- `PC5_BASELINE_STATUS := FROZEN_READY_FOR_FIRST_RUN`
- `L_B1_03_EXPORT_STATUS := NOT_EXPORTED`
- `B1_STATUS := OPEN`

## Decyzja profesorska (następny uczciwy krok)
Uruchomić `P1408_STRICT_ONLY_PC5_FIRST_DUAL_METRIC_RUN` na zamrożonych progach i opublikować jawny PASS/FAIL.

## Omówienie dla laika
To etap ustawienia nowej wersji modelu: najpierw zamrażamy reguły i progi, aby kolejny test był uczciwy i porównywalny.
