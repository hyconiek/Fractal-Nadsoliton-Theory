# RAPORT QW-1865: SYNTHETIC RECOVERY BENCHMARK

- Data UTC: 2026-03-03T00:04:01.329166+00:00
- Verdict: **SYNTHETIC_RECOVERY_PARTIAL**
- n_rep: 240

## Feature Set
- phase_increment, envelope_decay, zero_cross_offset

## Metrics
- RMSE omega/phi/beta: 0.0197 / 0.0384 / 0.0148
- median abs error omega/phi/beta: 0.0134 / 0.0265 / 0.0079
- tolerance hit omega(|e|<=0.08): 1.000
- tolerance hit phi(|e|<=0.20): 1.000
- tolerance hit beta(|e|<=0.015): 0.775
- nonboundary rate: 0.742
- error corr(omega,beta): -0.121

## Artifacts
- JSON: `report_qw1865_synthetic_recovery_benchmark.json`
