# RAPORT QW-1727: TOE EMPIRICAL CLOSURE GATE

- Data UTC: 2026-03-02T16:38:37.172662+00:00
- Global score: 0.143
- Readiness: **TOE_OPEN_NOT_EMPIRICALLY_CLOSED**
- Hard gate: **TOE_HARD_GATE_FAIL**

## Wyniki domenowe
- Flavor Extended Operator: FAIL (0.00) | FLAVOR_EXTENDED_OPERATOR_NOT_CLOSED
- Mass Strict OOS: PASS (1.00) | STRICT_OOS_PASS
- Parameter Stability: FAIL (0.00) | PARAMETERS_UNSTABLE_UNDER_PERTURBATION
- Unified Effective Hamiltonian: FAIL (0.00) | UNIFIED_EFFECTIVE_HAMILTONIAN_NOT_CLOSED
- GW Method Audit: FAIL (0.00) | risk_points=17
- GW Strict Reanalysis: FAIL (0.00) | GW_CROSS_HURST_ANOMALY_NOT_ROBUST
- GW FIN Projection: FAIL (0.00) | FIN_023_TO_031_PROJECTION_NOT_SUPPORTED

## Core gate flags
- Flavor Extended Operator: False
- Mass Strict OOS: True
- Parameter Stability: False
- Unified Effective Hamiltonian: False
- GW Strict Reanalysis: False
- GW FIN Projection: False

## Kolejne badania
- QW-1728: Latency-aware cross-detector analysis (explicit 10ms propagation model in estimator).
- QW-1729: Blind holdout epoch benchmark on untouched GPS windows.
- QW-1730: Unified symbolic derivation replacing fitted mapping lambda->flavor operator.
- QW-1731: External replication package (single-command reproducibility, frozen seeds, frozen manifests).

## Artefakty
- JSON szczegolowy: `report_qw1727_toe_empirical_closure_gate.json`
