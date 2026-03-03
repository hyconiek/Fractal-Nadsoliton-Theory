# RAPORT QW-1907: PRE-1700 TUNING BOUNDARY AUDIT (QW-700..1699)

## Wynik
- Analysis-parameter tuning pre-1700: **DETECTED**
- Kernel-core retuning pre-1700 (static signal): **NO_EXTERNAL_DATA_RETUNING_SIGNAL**
- Overall: **PRE1700_HAS_INFERENCE_BUT_NO_EXTERNAL_KERNEL_RETUNING_SIGNAL**

## Zakres i statystyki
- files_scanned: 405
- rows_total: 405
- rows_with_inference: 74
- rows_inference_simulation: 18
- rows_inference_external: 0
- rows_with_external_markers: 26
- rows_with_freeze_claim: 40
- rows_kernel_param_sweep_candidates: 0
- rows_kernel_param_sweep_candidates_with_external_data: 0

## Pierwsze wystapienia (QW id)
- inference_any: 713
- inference_simulation: 747
- inference_external: None
- external_any: 1500
- kernel_param_sweep_candidate: None

## Przykladowe pliki (head)
- inference_simulation:
  - QW_1525_Observational_Test_Pipeline.py
  - QW-747_to_QW-766_Batch_Suite.py
  - QW_1529_Rubikon_Test.py
  - QW_1526_Standalone_Simulation.py
  - QW_1528b_Relative_Test.py
  - QW_1533_Rubikon_Final_Audit.py
  - QW_1503_Dark_Energy_Pruning.py
  - QW-827_to_QW-846_True_Nadsoliton_Suite.py
  - QW_1557_Transport_Audit.py
  - QW_1501_Generative_Vacuum.py
- inference_external:
- kernel_param_sweep_candidates:

## Ograniczenie
- Static lexical audit. Detects markers, not full causal proof. Use as screening/audit evidence, not as final semantic adjudication.

- JSON: `report_qw1907_pre1700_tuning_boundary_audit.json`
