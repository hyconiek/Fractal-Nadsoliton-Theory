# P1411 — STRICT-ONLY PC6 PROVIDER BASELINE DESIGN (PL)

## Cel kroku
Po aktywacji `PC6` w `P1410` zamrażamy nową bazę provider-class dla strict-only toru
`F_Nadsoliton ⇒ L_SM + L_GR` bez bridge do legacy.

## Zakres strict-only
- `legacy_bridge_used = false`
- `inherits_from_pc5_loop = false`
- `noncyclic_provider_shift = enforced`

## Baseline PC6 (v1)
- `provider_class_id = PC6_strict_phase_bridge_selector_v1`
- `anchor_family = A6_noncyclic_phase_bridge_stabilizers`
- `epsilon_sign_v1 = 0.05`
- `epsilon_drift_v1 = 0.04`
- `first_run_contract = dual_metric_required`

## Status kroku
- `PC6_BASELINE_STATUS := FROZEN_READY_FOR_FIRST_RUN`
- `L_B1_03_EXPORT_STATUS := NOT_EXPORTED`
- `B1_STATUS := OPEN`

## Decyzja profesorska (następny uczciwy krok)
Uruchomić `P1412_STRICT_ONLY_PC6_FIRST_DUAL_METRIC_RUN` i opublikować jawny PASS/FAIL bez retuningu progów.

## Omówienie dla laika
To przygotowanie nowej wersji modelu: ustawiamy stabilne reguły testu, żeby wynik kolejnego kroku był rzetelny i porównywalny.
