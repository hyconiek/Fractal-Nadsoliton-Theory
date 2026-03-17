# P704 Current Non-Strict: Standard Model Host Matching From `F704` H\_psi Eigenvalue Proxy (No False-PASS)

Status: `P704_CURRENT_NONSTRICT_STANDARD_MODEL_HOST_MATCHING_FROM_F704_H_PSI_EIGENVALUE_PROXY_PROBE_NO_FALSE_PASS`  
As of: `2026-03-17`

## Goal

Compute **non-strict** host-matching metrics between:

1. the strict **basis-invariant** mass proxy spectrum exported by `F704` (eigenvalues of `H_psi`), and
2. an explicit external Standard Model mass target dataset (`external_data/sm_mass_targets_v1.json`),
   under an explicit **versioned** non-strict host-matching policy (`external_data/sm_host_matching_policy_v1.json`).

This probe is **outside strict scope** by construction: it depends on an external dataset and an identification/matching policy.  
It does not claim any Standard Model match or ToE closure.

## Inputs

- `F704` strict mass observable:
  - `fundamental_action_reconstruction/generated/mass_observable_diagonal_local_strict_derived_v1.json`
- External dataset (targets only):
  - `fundamental_action_reconstruction/external_data/sm_mass_targets_v1.json`
- External policy (identification / assignment / scale DOF):
  - `fundamental_action_reconstruction/external_data/sm_host_matching_policy_v1.json`

## Output

Executed by:

```bash
python3 fundamental_action_reconstruction/p704_current_nonstrict_standard_model_host_matching_from_f704_h_psi_eigenvalue_proxy_probe.py
```

Exports:

- `fundamental_action_reconstruction/generated/p704_current_nonstrict_standard_model_host_matching_from_f704_h_psi_eigenvalue_proxy_probe.json`
- `fundamental_action_reconstruction/generated/p704_current_nonstrict_standard_model_host_matching_from_f704_h_psi_eigenvalue_proxy_probe_summary.json`

## Hard limits

- Non-strict only (external dataset + matching policy).
- No physical unit calibration claim.
- No Standard Model identification claim.
- No `QW-2191` discharge.
- No ToE closure.
