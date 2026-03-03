# RAPORT QW-1866: PAIRED SIGNED-INTERVENTION BENCHMARK

- Data UTC: 2026-03-03T00:04:01.866073+00:00
- Verdict: **PAIRED_SIGNED_INTERVENTION_NOT_SUPPORTED**
- n_rep: 180

## Feature Set
- phase_increment, envelope_decay, zero_cross_offset

## Baseline vs Paired RMSE
- omega: 0.0199 -> 0.0206
- phi: 0.0386 -> 0.0267
- beta: 0.0144 -> 0.0212

## Baseline vs Paired Tolerance Hit
- omega(|e|<=0.08): 1.000 -> 1.000
- phi(|e|<=0.20): 1.000 -> 1.000
- beta(|e|<=0.015): 0.733 -> 0.639

## Coupling Diagnostic
- corr(omega_error, beta_error): -0.114 -> -0.156

## Improvement Factors
- rmse omega factor: 0.966
- rmse phi factor: 1.445
- rmse beta factor: 0.679
- beta tolerance gain: -0.094
- |corr| reduction: -0.042

## Artifacts
- JSON: `report_qw1866_paired_signed_intervention_benchmark.json`
