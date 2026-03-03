# RAPORT QW-1731: NODE COMPATIBILITY AUDIT

- Data UTC: 2026-03-02T16:43:36.134939+00:00
- Verdict: **NODE_NARRATIVE_STRONGLY_INCONSISTENT**

## Priors from Nadsoliton Characteristics
- alpha_geo = 2.772588722240
- omega = 0.785398163397 (pi/4)
- phi = 0.523598775598 (pi/6)
- beta_tors = 0.01000

## Analytic Zero Sequence from Priors
- first_zero = 1.3333
- zero_spacing = 4.0000

## Claim Set A: nodes [2,5,8,11]
- prior_loss = 0.500000
- best_fit omega = 1.046857, phi = -0.521504, loss = 1.317610e-06
- relative shift from priors: domega=0.333, dphi/pi=0.333
- mean distance to prior zero sequence = 1.0000

## Claim Set B: nodes [2,8,14]
- prior_loss = 0.416667
- best_fit omega = 0.523571, phi = 0.524646, loss = 7.042304e-07
- relative shift from priors: domega=0.333, dphi/pi=0.000
- mean distance to prior zero sequence = 0.8889

## Hierarchy Check
- abs(K(7))/abs(K(1)) with priors = 3.523

## Flags
- NODE_SET_A_INCOMPATIBLE_WITH_CHARACTERISTIC_PRIORS
- NODE_SET_B_INCOMPATIBLE_WITH_CHARACTERISTIC_PRIORS
- NODE_SET_A_FAR_FROM_ANALYTIC_ZERO_SEQUENCE
- NODE_SET_B_FAR_FROM_ANALYTIC_ZERO_SEQUENCE

## Artifacts
- JSON: `report_qw1731_nadsoliton_kernel_node_compatibility.json`
