# N693 Current Strict `T173` Output Sign‑Lift Gauge Covariance Theorem

Status: `N693_CURRENT_STRICT_T173_OUTPUT_SIGN_LIFT_GAUGE_COVARIANCE_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-17`

## Goal

Package one theorem-level statement that the directed closure output sign‑lift is a tracked `Z2` gauge datum and transforms covariantly under chartwise sign relifts of the directed state representative.

## Inputs

1. Premise-based directed closure object on `C_v1` with explicit output sign‑lift:
   - `SelectorClosure_global_C_v1_directed_strict_v1` (`F677`).
2. Sign-fixed directed state representative with explicit chartwise sign relift `t_by_pair`:
   - `SelectorState_global_C_v1_directed_sign_fixed_from_strict_core_payload_weights_strict_convention_v1` (`F690`).
3. Sign-fixed directed closure object on `C_v1` with explicit output sign‑lift:
   - `SelectorClosure_global_C_v1_directed_sign_fixed_from_strict_core_payload_weights_strict_convention_v1` (`F692`).
4. Gauge-covariance audit:
   - `P693` audits the identity `s_out^fix(pair) = t(pair) · s_out^prem(pair)` on `{pair1..pair5}`.

## Theorem (scope-limited)

On the current repo state, the directed closure output sign‑lift satisfies chartwise `Z2` gauge covariance:

```text
s_out^fix(pair_m) = t(pair_m) * s_out^prem(pair_m)   for m=1..5,
```

where `t` is the exported chartwise sign relift mapping the premise-based directed representative to the exported sign-fixed representative.

Therefore, output sign‑lift data is convention/gauge and does not upgrade to any strict physical sign datum.

## Hard limits

This theorem does **not** claim:

- any directed/sign-sensitive physical orientation datum in strict core,
- kernel-alone/global `QW-2191` discharge,
- ToE closure.

## Exported artifact

- `fundamental_action_reconstruction/generated/n693_current_strict_t173_output_sign_lift_gauge_covariance_theorem_summary.json`

