# P478 Current Strict `T169` `r_ordpow` Sign‑Scan for `R18` Pair1 Zero‑Equations Probe

Status: `P478_EXECUTED_T169_RORDPOW_SIGN_SCAN_FOR_R18_PAIR1_ZERO_EQUATIONS_PROBE_NO_FALSE_PASS`  
As of: `2026-03-16`

## Goal

After `R18`, the missing diagonal-side ingredient for the existing-kernel-feedback host matching route is a strict
**zero/cancellation witness** for the declared `pair1` residual block.

`P477/N520` already show that the currently exported strict-derived value instance (from `P459`, using conditional
`N477`) violates the `R18` declared `pair1` residual zero equations (`c1c1=0`, `c1s1=0`, `s1s1=0`).

This probe answers one narrower “no false-pass” question:

```text
Is the violation merely an artifact of the specific T169 sign-selection choice in F447,
or does it persist for all sign vectors under the same fixed T169 r_ordpow magnitude lift,
when evaluated under the same conditional N477 rewrite?
```

So `P478` brute-forces the finite sign space:

- scan all `2^11=2048` sign vectors (fixing `s0=+1`),
- keep magnitudes fixed to the strict `r_ordpow` magnitude lift,
- keep `g4` fixed to the uniform T169 lift,
- compute `d_k` via the conditional `N477` Yukawa-free diagonal residual rewrite,
- evaluate the `R18` declared `pair1` zero equations.

## Inputs

1. `R14` strict frozen kernel specialization (`K_total`):
   - `fundamental_action_reconstruction/generated/r14_explicit_frozen_kernel_specialization_packet_for_host_matching_route.json`
2. `R15` diagonal floor `m0^2`:
   - `fundamental_action_reconstruction/generated/r15_explicit_host_scalar_floor_embedding_packet_for_host_matching_route.json`
3. `R18` pair1 coefficient-class reduction and exact zero system:
   - `fundamental_action_reconstruction/generated/r18_pair1_residual_declared_pullback_coefficient_class_reduction_packet_for_host_matching_route.json`
4. `F447` strict-derived T169 provider packet (as the source of the fixed `r_ordpow` magnitudes and scalar constants):
   - `fundamental_action_reconstruction/generated/f447_current_strict_t169_qw2122_scalar_to_t168_per_site_value_provider_strict_derived_v1.json`

## Output artifacts

- `fundamental_action_reconstruction/generated/pair1_residual_zero_equations_sign_scan_under_rordpow_magnitudes_value_instance_v1.json`
- `fundamental_action_reconstruction/generated/p478_current_strict_t169_rordpow_sign_scan_for_r18_pair1_zero_equations_probe.json`
- `fundamental_action_reconstruction/generated/p478_current_strict_t169_rordpow_sign_scan_for_r18_pair1_zero_equations_probe_summary.json`

## Hard limits (no false pass)

This probe does **not** claim:

1. any strict vacuum stationarity witness (it uses conditional `N477`),
2. any strict zero witness for the declared `pair1` residual equations (unless later promoted by a theorem-level discharge),
3. host matching (`C10_B1`) discharge,
4. selector-relevant physical canonicalization inside the `QW-2191` `O(2)` family,
5. strict-core selector closure / admissible `S_sel_int`,
6. global `QW-2191` discharge,
7. ToE closure.

