# RAPORT QW-1840: JOINT CONFIRMATORY POWER SENSITIVITY

- Data UTC: 2026-03-02T22:34:58.638666+00:00
- Verdict: **JOINT_CONFIRMATORY_HIGH_POWER**
- Protocol SHA256 (QW-1839): `b9cf21d3d32508e95c6f7cef2e8a953a12f6c7ea7732b6288605c518ab5db5af`
- Rozmiary bootstrap: PTA n=14, GW n=5

## Marginesy do progow (stan referencyjny)
- PTA mean margin: 0.022399
- PTA std margin: 0.009273
- GW AUC margin: 0.041492
- GW ADV margin: 0.157098
- GW GAP margin: 0.000346

## Wyniki scenariuszy (consumed margin)
- margin_consumption_0pct: p(PTA)=0.994 [0.993,0.995], p(GW)=1.000 [1.000,1.000], p(JOINT)=0.994 [0.993,0.995]
- margin_consumption_50pct: p(PTA)=0.698 [0.693,0.703], p(GW)=0.998 [0.998,0.999], p(JOINT)=0.696 [0.691,0.702]
- margin_consumption_80pct: p(PTA)=0.505 [0.499,0.510], p(GW)=0.843 [0.839,0.847], p(JOINT)=0.427 [0.421,0.432]
- margin_consumption_100pct: p(PTA)=0.327 [0.321,0.332], p(GW)=0.294 [0.289,0.299], p(JOINT)=0.095 [0.092,0.098]
- margin_consumption_120pct: p(PTA)=0.157 [0.153,0.161], p(GW)=0.027 [0.025,0.029], p(JOINT)=0.004 [0.003,0.005]

## Robustnosc joint gate
- Max consumed margin z p(JOINT)>=0.90: 0.20
- Max consumed margin z p(JOINT)>=0.50: 0.75

## Artifacts
- JSON: `report_qw1840_joint_confirmatory_power_sensitivity.json`
