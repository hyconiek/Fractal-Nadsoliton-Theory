# RAPORT QW-1718: PROSPECTIVE OOS BENCHMARK

- Data UTC: 2026-03-02T15:45:32.390370+00:00
- Manifest hash: `d5815383290f5d387bd5baa744a94b770a27054fd8153d149b08f5e350b80c1b`
- Score: 0.600 (3/5)
- Werdykt: **PROSPECTIVE_BENCHMARK_PARTIAL**

## Checki
- Flavor shared CKM/PMNS: pass=False, value={'ckm': 21.27155345397912, 'pmns': 33.670736737640965}, threshold=15.0
- Mass OOS test_mean: pass=True, value=8.999014915803762, threshold=10.0
- Mass bootstrap robustness: pass=False, value={'test_median_pct': 10.830672556718213, 'gap_mean_pct': 21.20328180351618}, threshold={'test_median_pct': 10.0, 'gap_mean_pct': 7.0}
- IR Lorentz delta_max: pass=True, value=1.147288708480687e-09, threshold=0.001
- Mass identifiability cond(X^T X): pass=True, value=4.761438831454202, threshold=100000.0

## Artefakty
- JSON szczegółowy: `report_qw1718_prospective_oos_benchmark.json`
