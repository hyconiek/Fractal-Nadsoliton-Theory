# P403 Current Strict Sigma-Int Theta-Pair Delta_d Sensitivity Audit Probe

Status: `P403_EXECUTED_CURRENT_STRICT_SIGMA_INT_THETA_PAIR_DELTA_D_SENSITIVITY_AUDIT_PROBE_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

`F325/N436` export one strict-side **candidate** theta-pair selector ingredient:

```text
ThetaPair_sigma_int_strict_selector_ingredient_o2_cut_candidate_v1
```

However, the underlying positive-window construction (T119-style) still
contains one explicit free slot:

```text
delta_d ∈ (0, d_local/11]
```

and different admissible `delta_d` choices produce different theta-pair
outputs.

This probe quantifies that dependence on the current strict tuple and strict
sigma-int provenance, without promoting any result to:

- strict-core theta export,
- actual object-support above the export-map object (`N302/N395`),
- admissible `S_sel_int` / selector closure,
- `QW-2191` discharge,
- ToE closure.

## Implementation

Executed by:

```text
python3 fundamental_action_reconstruction/p403_current_strict_sigma_int_theta_pair_delta_d_sensitivity_audit_probe.py
```

Outputs:

```text
fundamental_action_reconstruction/generated/p403_current_strict_sigma_int_theta_pair_delta_d_sensitivity_audit_probe.json
fundamental_action_reconstruction/generated/p403_current_strict_sigma_int_theta_pair_delta_d_sensitivity_audit_probe_summary.json
```

## Exact verdict

On the current repo state:

```text
P403: PASS (audit executed).
The theta-pair output varies over admissible delta_d choices inside the positive-window corridor.
Therefore, without an explicit delta_d selector premise, the construction remains candidate-only
and cannot be promoted to strict-core theta export or to actual object-support discharge.
```

No false pass: this probe is a sensitivity audit only.

