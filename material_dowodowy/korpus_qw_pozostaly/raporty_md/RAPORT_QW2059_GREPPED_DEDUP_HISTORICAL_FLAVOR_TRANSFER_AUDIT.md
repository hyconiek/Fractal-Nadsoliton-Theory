# RAPORT QW-2059: GREPPED DEDUP + HISTORICAL FLAVOR TRANSFER AUDIT

- Data UTC: 2026-03-04T02:04:11.971828+00:00
- Verdict: **DEDUP_AUDIT_IDENTIFIES_EXISTING_METHODS_AND_NO_STRICT_PASS_UNDER_CURRENT_KERNEL**

## Grep Dedup Index
- indexed files count: 57
- QW_1966_isospin_split_scan: True
- QW_2029_shared_flavor_basis_scan: True
- QW_2012_no_ansatz_no_fit: True
- QW_2056_operator_family_frontier: True
- QW_2057_su3_rotation_frontier: True
- QW_2058_nonabelian_no_fit: True

## Evaluations
### transfer_qw1966_best_to_qw2049
- method_family: isospin_split_shared_flavor_dynamics_scan
- flavor CKM/PMNS mean rel%: 41.312/14.780
- GW auc/adv/sep/gap: 0.8421/0.4012/0.003790/0.003101
- pass_count: 4/7
- all_pass: False
- notes: historical scan-derived parameters, not no-fit first-principles

### transfer_qw2029_best_to_qw2049
- method_family: shared_flavor_basis_scan
- flavor CKM/PMNS mean rel%: 11.867/9.386
- GW auc/adv/sep/gap: 0.8427/0.4012/0.003868/0.003124
- pass_count: 5/7
- all_pass: False
- notes: historical scan-derived parameters, not no-fit first-principles

### reference_qw2058_nonabelian_no_fit
- method_family: nonabelian_no_fit
- flavor CKM/PMNS mean rel%: 61.378/68.638
- GW auc/adv/sep/gap: 0.8320/0.3419/0.003108/0.002611
- pass_count: 4/7
- all_pass: False
- notes: current strict first-principles baseline

## Required Next Step
- BUILD_QW2060_FIRST_PRINCIPLES_RECONSTRUCTION_OF_QW2029_FLAVOR_BASIS_FROM_KERNEL_INVARIANTS_ONLY

## Artifacts
- JSON: `report_qw2059_grepped_dedup_historical_flavor_transfer_audit.json`
