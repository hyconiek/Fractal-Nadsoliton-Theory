# RAPORT QW-2056: FIRST-PRINCIPLES FLAVOR OPERATOR FAMILY FRONTIER

- Data UTC: 2026-03-04T01:49:37.709710+00:00
- Verdict: **FIRST_PRINCIPLES_OPERATOR_FAMILY_FRONTIER_FAILS_TO_CLOSE_FLAVOR**
- any_all_pass: False

## Best Row
- family: phase_sign
- q_scheme: proxy_old
- pass_count: 3/6
- CKM fixed / PMNS / CKM perm-envelope: 62.308 / 32.138 / 48.538
- GW auc/adv/sep/gap: 0.8485/0.4130/0.004815/0.003727

## Best Flavor Fixed (CKM+PMNS)
- family/q_scheme: phase_sign / proxy_old | CKM=62.308 | PMNS=32.138

## Best CKM (Fixed Order)
- family/q_scheme: locality / proxy_old | CKM=57.084 | PMNS=43.937

## Best CKM (Permutation Envelope, Diagnostic)
- family/q_scheme: critical / proxy_old | CKM_perm_env=24.337 | PMNS=36.400

## Top 10 Rows
- 1. phase_sign | proxy_old | pass 3/6 | CKM=62.308 | PMNS=32.138 | CKM_perm_env=48.538
- 2. legacy | proxy_old | pass 4/6 | CKM=197.086 | PMNS=14.060 | CKM_perm_env=66.041
- 3. critical | proxy_old | pass 3/6 | CKM=104.919 | PMNS=36.400 | CKM_perm_env=24.337
- 4. locality | proxy_old | pass 3/6 | CKM=57.084 | PMNS=43.937 | CKM_perm_env=46.630
- 5. ultra_local | proxy_old | pass 3/6 | CKM=63.398 | PMNS=45.720 | CKM_perm_env=50.152
- 6. phase_sign | quark_mass_inversion | pass 3/6 | CKM=602.034 | PMNS=32.138 | CKM_perm_env=53.778
- 7. ultra_local | quark_mass_inversion | pass 3/6 | CKM=625.386 | PMNS=45.720 | CKM_perm_env=42.149
- 8. locality | quark_mass_inversion | pass 3/6 | CKM=656.923 | PMNS=43.937 | CKM_perm_env=32.813
- 9. critical | quark_mass_inversion | pass 3/6 | CKM=829.877 | PMNS=36.400 | CKM_perm_env=103.817
- 10. legacy | quark_mass_inversion | pass 4/6 | CKM=987.114 | PMNS=14.060 | CKM_perm_env=189.005

## Required Next Step
- DERIVE_NONABELIAN_FLAVOR_GENERATOR_FROM_KERNEL_DYNAMICS_WITHOUT_NEW_FREE_PARAMETERS

## Artifacts
- JSON: `report_qw2056_first_principles_flavor_operator_family_frontier.json`
