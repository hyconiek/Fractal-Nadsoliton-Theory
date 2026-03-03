# RAPORT QW-1881: NODE-STATE IDENTIFIABILITY STRESS

- Data UTC: 2026-03-03T00:35:37.484438+00:00
- Verdict: **NODE_STATE_IDENTIFIABILITY_STRESS_FAIL**
- source profiles: 6
- stress samples: 1080

## Global Metrics
- omega_iqr: 0.0808
- phi_circular_std: 0.0852
- beta_iqr: 0.0150
- corr(omega,beta): -0.2434
- nonboundary_rate: 0.126
- stability_index: 0.544

## Scenario Summary
- noise_low: n=216, omega_iqr=0.0565, beta_iqr=0.0157, canon_median=0.9178, nonboundary=0.079
- noise_mid: n=216, omega_iqr=0.0944, beta_iqr=0.0149, canon_median=0.9267, nonboundary=0.111
- noise_high: n=216, omega_iqr=0.1115, beta_iqr=0.0172, canon_median=0.9362, nonboundary=0.171
- edge_mid: n=216, omega_iqr=0.0601, beta_iqr=0.0163, canon_median=0.9170, nonboundary=0.093
- combined_high: n=216, omega_iqr=0.0779, beta_iqr=0.0127, canon_median=0.9418, nonboundary=0.176

## Artifacts
- JSON: `report_qw1881_node_state_identifiability_stress.json`
