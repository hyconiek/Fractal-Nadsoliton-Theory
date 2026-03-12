# P407 Current Strict Sigma-Int Theta-Pair Eps Sensitivity Audit Probe

Status: `P407_EXECUTED_CURRENT_STRICT_SIGMA_INT_THETA_PAIR_EPS_SENSITIVITY_AUDIT_PROBE_NO_FALSE_PASS`  
As of: `2026-03-12`

## Goal

`F325/N436` export one strict-side **candidate** theta-pair selector ingredient:

```text
ThetaPair_sigma_int_strict_selector_ingredient_o2_cut_candidate_v1
```

However, the underlying sigma-int → `E_pair` generator contract (`T117`) carries
one explicit free slot:

```text
eps ∈ [0,1]
```

and different admissible `eps` choices produce different theta-pair outputs.

This probe quantifies that dependence on the current strict tuple and strict
sigma-int provenance, at fixed corridor-saturation `delta_d`, without promoting
any result to:

- strict-core theta export,
- actual object-support above the export-map object (`N302/N395`),
- admissible `S_sel_int` / selector closure,
- `QW-2191` discharge,
- ToE closure.

## Implementation

Executed by:

```text
python3 fundamental_action_reconstruction/p407_current_strict_sigma_int_theta_pair_eps_sensitivity_audit_probe.py
```

Outputs:

```text
fundamental_action_reconstruction/generated/p407_current_strict_sigma_int_theta_pair_eps_sensitivity_audit_probe.json
fundamental_action_reconstruction/generated/p407_current_strict_sigma_int_theta_pair_eps_sensitivity_audit_probe_summary.json
```

## Exact verdict

On the current repo state:

```text
P407: PASS (audit executed).
The theta-pair output varies over admissible eps choices in [0,1] (at fixed delta_d inside the corridor).
Therefore, without an explicit eps selector premise (or a strict-derived eps law), the construction remains
candidate-only and cannot be promoted to strict-core theta export or to actual object-support discharge.
```

No false pass: this probe is a sensitivity audit only.

