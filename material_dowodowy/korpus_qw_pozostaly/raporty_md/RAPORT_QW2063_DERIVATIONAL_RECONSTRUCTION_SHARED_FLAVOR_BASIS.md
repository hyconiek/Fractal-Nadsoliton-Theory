# RAPORT QW-2063: DERIVATIONAL RECONSTRUCTION SHARED FLAVOR BASIS

- Data UTC: 2026-03-04T02:19:38.612949+00:00
- Verdict: **DERIVATIONAL_RECONSTRUCTION_TRIAD_PASS_PHYSICAL_PROVISIONAL_FIRST_PRINCIPLES**
- pass_count: 11/12
- physical_pass: True

## Metrics
- mass mean/max/tau-charm rel%: 12.051/34.013/14.013
- flavor CKM/PMNS mean rel%: 11.867/9.386
- GW auc/adv/sep/gap: 0.8150/0.3103/0.002056/0.001289

## Derivation Snapshot
- invariants: {'decay_ratio': 0.10006511151920049, 'phi_abs': 0.1624999999999998, 'omega': 0.18575000000000005, 'eta': 1.8}
- q_nu_order: [2, 1, 0]
- params: {'p_amp': 0.7, 'r_dist': -0.2, 'lambda_mix': 0.8, 'rho_gap': 0.4, 'chi_im': 0.3, 'phase_q': 0.25, 'phase_q3': 0.0, 'theta_iso': 0.6, 'theta_sector': 0.3, 'diag_q_coeff': 0.1, 'amp_qbias': 0.4, 'diag_iso': 0.0, 'diag_sector': 0.0}

## Robustness
- pass_rate: 0.897 (269/300)

## Flags
- mass_mean_rel_pct_le_max: True
- mass_max_rel_pct_le_max: True
- mass_tau_charm_ratio_err_le_max: True
- ckm_mean_rel_pct_le_max: True
- pmns_mean_rel_pct_le_max: True
- gw_sep_ge_min: True
- gw_adv_ge_min: True
- gw_auc_ge_min: True
- gw_control_gap_le_max: True
- deterministic_no_scan: True
- local_robustness_pass_rate_ge_0p80: True
- strict_first_principles_foundational_constants_derived: False

## Required Next Step
- FORMALLY_DERIVE_REGIME_CONSTANTS_FROM_NADSOLITON_THEORY_TO_PROMOTE_TO_STRICT_FIRST_PRINCIPLES

## Artifacts
- JSON: `report_qw2063_derivational_reconstruction_shared_flavor_basis.json`
