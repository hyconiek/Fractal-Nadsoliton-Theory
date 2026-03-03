# RAPORT QW-1712: IR LORENTZ RECOVERY TEST

- Data UTC: 2026-03-02T15:33:21.114499+00:00
- Werdykt: **IR_LORENTZ_RECOVERY_OK**

## 1) Baseline (beta_tors = 0.01)
- max Delta_Lorentz (mu <= 1e-3): 1.147289e-09
- mean Delta_Lorentz (mu <= 1e-3): 8.357683e-11
- Delta_Lorentz (mu=1): 1.189196e-03

## 2) Sensitivity beta_tors
- beta=0.005: delta_ir_max=5.736442e-10, delta_ir_mean=4.178839e-11, delta_uv(mu=1)=5.861152e-04
- beta=0.010: delta_ir_max=1.147289e-09, delta_ir_mean=8.357683e-11, delta_uv(mu=1)=1.189196e-03
- beta=0.020: delta_ir_max=2.294577e-09, delta_ir_mean=1.671537e-10, delta_uv(mu=1)=2.446128e-03
- beta=0.050: delta_ir_max=5.736444e-09, delta_ir_mean=4.178844e-10, delta_uv(mu=1)=6.620743e-03

## 3) Interpretacja
- Test mówi, czy model EFT zachowuje zgodność Lorentzowską w limicie IR.
- Nie jest to pełny dowód relatywistycznej kompletności; to test konieczny, nie wystarczający.

## Artefakty
- JSON szczegółowy: `report_qw1712_ir_lorentz_recovery_test.json`
