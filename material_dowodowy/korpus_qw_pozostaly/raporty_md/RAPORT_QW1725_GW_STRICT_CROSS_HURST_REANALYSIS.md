# RAPORT QW-1725: GW STRICT CROSS-HURST REANALYSIS

- Data UTC: 2026-03-02T16:38:20.145913+00:00
- Werdykt: **GW_CROSS_HURST_ANOMALY_NOT_ROBUST**

## Baseline
- H_cross(H1,L1) @ N=524288: 0.3143
- H_cross(H1,V1) @ N=524288: 0.7603605018723744

## Null models
- Phase-randomized: mean=0.2851, std=0.0114, p_lower=1, z=2.57
- Circular-shift: mean=0.3045, std=0.0110, p_lower=0.8017, z=0.89

## Stability
- length spread: 0.0990
- window std: 0.0126

## Lag test (+/-25 ms)
- corr(+10 ms): 0.00163
- best abs lag: -4.39 ms (corr=0.00262)

## Artefakty
- JSON szczegolowy: `report_qw1725_gw_strict_cross_hurst_reanalysis.json`
