# P693 Current Strict `T173` Output Sign‑Lift Gauge Covariance Audit Probe (No False‑PASS)

Status: `P693_CURRENT_STRICT_T173_OUTPUT_SIGN_LIFT_GAUGE_COVARIANCE_AUDIT_PROBE_NO_FALSE_PASS`  
As of: `2026-03-17`

## Goal

Audit a precise gauge‑covariance statement for directed closure outputs:

Let:

- `s_out^prem(pair_m)` be the per‑chart output sign‑lift used by the premise‑based directed closure object
  `SelectorClosure_global_C_v1_directed_strict_v1` (`F677`), and
- `t(pair_m)` be the per‑chart sign relift mapping the premise‑based directed state representative to the exported sign‑fixed directed state representative (`F690`),
- `s_out^fix(pair_m)` be the per‑chart output sign‑lift used by the sign‑fixed directed closure object (`F692`).

Then, because `v_out_raw` transforms by the same chartwise sign relift under `u -> t·u`, we must have:

```text
s_out^fix(pair_m) = t(pair_m) * s_out^prem(pair_m)   for m=1..5.
```

This probe audits the identity on exported artifacts, confirming that output sign‑lift data behaves as tracked `Z2` gauge, not as a strict physical sign datum.

## Hard limits

This probe does **not** claim:

- any directed/sign‑sensitive physical orientation datum in strict core,
- kernel-alone/global `QW-2191` discharge,
- ToE closure.

