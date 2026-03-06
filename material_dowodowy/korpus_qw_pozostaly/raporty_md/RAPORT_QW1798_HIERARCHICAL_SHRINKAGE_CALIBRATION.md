# RAPORT QW-1798: HIERARCHICAL SHRINKAGE CALIBRATION

- Data UTC: 2026-03-02T20:47:34.084914+00:00
- Family: mixed_low [1, 2, 3]
- Best shrink factor: 0.20 (CONDITIONAL)
- Full delta hier-M2 (best): 2.0490
- P(hier>M2) (best): 1.000
- P(hier>flat) (best): 1.000
- Std delta hier-M2 (best): 0.422
- Residual improvements (best): d|rho|=0.1500, d(bin)=-0.1053
- Verdict: **HIERARCHICAL_SHRINKAGE_PARTIAL**

## Sweep
- sf=0.20: full_delta=2.0490, P(hier>M2)=1.000, P(hier>flat)=1.000, std_delta=0.422, score=0.800, pass_basic=False
- sf=0.35: full_delta=2.0529, P(hier>M2)=1.000, P(hier>flat)=1.000, std_delta=0.576, score=0.800, pass_basic=False
- sf=0.50: full_delta=1.9876, P(hier>M2)=1.000, P(hier>flat)=1.000, std_delta=0.457, score=0.800, pass_basic=False
- sf=0.70: full_delta=1.8954, P(hier>M2)=1.000, P(hier>flat)=1.000, std_delta=0.475, score=0.800, pass_basic=False
- sf=1.00: full_delta=1.7164, P(hier>M2)=1.000, P(hier>flat)=1.000, std_delta=0.589, score=0.800, pass_basic=False

## Artifacts
- JSON: `report_qw1798_hierarchical_shrinkage_calibration.json`
