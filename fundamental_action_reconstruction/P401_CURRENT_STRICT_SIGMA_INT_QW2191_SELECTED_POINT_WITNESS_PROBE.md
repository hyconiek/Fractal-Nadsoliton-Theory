# P401 Current Strict Sigma-Int QW-2191 Selected-Point Witness Probe

Status: `P401_EXECUTED_CURRENT_STRICT_SIGMA_INT_QW2191_SELECTED_POINT_WITNESS_PROBE_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

After `F325/N436`, the repo exports one strict-side **candidate** selector-ingredient object:

```text
ThetaPair_sigma_int_strict_selector_ingredient_o2_cut_candidate_v1
```

persisted as:

```text
fundamental_action_reconstruction/generated/theta_pair_sigma_int_strict_selector_ingredient_o2_cut_candidate_v1.json
```

This probe adds one missing sanity link to the `QW-2190/QW-2191` scaffold:

```text
do the selected candidate angles (theta_1^cand, theta_2^cand)
actually lie on the QW-2191 exhibited O(2) rotation family,
i.e. do they preserve the same kernel-subspace invariance and Lie-closure audits?
```

This is a **selected-point witness only** (compatibility check), not a uniqueness proof.

## Strict-admissible evidence reused

1. `F325/N436`
   - exports the candidate theta-pair artifact (declared O(2)-cut ingredient; no closure claim).
2. `QW-2190`
   - fixed strict kernel tuple and deterministic real Fourier scaffold `(e0,c1,s1,c2,s2)` on `n=12`.
3. `QW-2191`
   - defines the kernel-subspace invariance residual and SU(3)/SU(2) Lie-closure residual audits and proves
     the continuous O(2) nonuniqueness obstruction from kernel algebra alone.

## Probe implementation artifact

The probe is executed by:

```text
python3 fundamental_action_reconstruction/p401_current_strict_sigma_int_qw2191_selected_point_witness_probe.py
```

and persists its results as:

```text
fundamental_action_reconstruction/generated/p401_current_strict_sigma_int_qw2191_selected_point_witness_probe.json
fundamental_action_reconstruction/generated/p401_current_strict_sigma_int_qw2191_selected_point_witness_probe_summary.json
```

## Exact verdict

On the current repo state:

```text
P401: PASS (selected-point compatibility witness).
The selected point (theta_1^cand,theta_2^cand) satisfies the same QW-2191 invariance/closure audits,
so it lies on the exhibited O(2) family as a declared representative point.
```

## What P401 proves

`P401` proves only:

1. the exported candidate theta pair is numerically compatible with the `QW-2191` O(2)-family audits
   (kernel invariance + Lie closure) at the selected point,
2. the exported `u_i^cand` vectors match the corresponding rotated Fourier basis vectors at the exported `theta_i^cand`
   (internal consistency of the packaging).

## What P401 does not prove

`P401` does not prove:

1. discharge of `QW-2191` (kernel-alone uniqueness obstruction remains true),
2. strict-core selector closure (no admissible strict `S_sel_int` is exported),
3. physical uniqueness of mode-index assignment from strict-core sources,
4. ToE closure.

No false pass: this probe is compatibility-only, not a uniqueness upgrade.

