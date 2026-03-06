# RAPORT QW-1732: FRACTAL PATH -> HYPERBOLIC TEST

- Data UTC: 2026-03-02T16:42:51.651808+00:00
- Total models scanned: 11011
- Best fit: nu=0.033, lam=0.067, beta=0.00424
- Best RMSE=1.960958e-03, R2=0.977902
- Coverage R2>=0.995: 0.0115
- Coverage R2>=0.98: 0.0326
- Verdict: **HYPERBOLIC_REDUCTION_PLAUSIBLE_BUT_TUNED**

## Ridge Statistics (R2>=0.995)
- count: 127
- mean(lam-nu): 0.6412073490813646
- std(lam-nu): 0.1362192081422473
- mean(beta): 0.13988993131177752
- std(beta): 0.04459417937924048

## Interpretation
- If high-quality region is narrow, path->hyperbolic claim is mechanistically plausible but not uniquely derived.
- If high-quality region is broad, reduction is robust under model perturbations.

## Artifacts
- JSON: `report_qw1732_fractal_path_to_hyperbolic_test.json`
