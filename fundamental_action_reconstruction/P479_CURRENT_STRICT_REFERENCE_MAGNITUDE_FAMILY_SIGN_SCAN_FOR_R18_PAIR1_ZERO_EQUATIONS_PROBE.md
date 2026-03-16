# P479 Current Strict Reference Magnitude Family Sign‑Scan for `R18` Pair1 Zero‑Equations Probe

Status: `P479_EXECUTED_REFERENCE_MAGNITUDE_FAMILY_SIGN_SCAN_FOR_R18_PAIR1_ZERO_EQUATIONS_PROBE_NO_FALSE_PASS`  
As of: `2026-03-16`

## Goal

After `R16–R18`, the diagonal-side bottleneck for the existing-kernel-feedback host matching / legacy chart-reduced operator
route (`P16`) is reduced to proving the **declared** `pair1` residual zero equations (`c1c1=0`, `c1s1=0`, `s1s1=0`),
plus the separate `QW-2191` selector-relevant canonicalization boundary.

`P477/N520` record that the currently exported strict-derived value instance (via `P459`, consuming the conditional `N477`
rewrite) violates the declared `R18` `pair1` zero equations.

`P478/N521` strengthen this: under the fixed strict `T169` `r_ordpow` magnitude lift (and fixed uniform `g4` lift),
no sign choice yields a solution to all three declared `R18` `pair1` zero equations under the same conditional `N477`
diagonal residual rewrite.

This probe answers the next honest narrow question:

```text
If we keep the same conditional N477 rewrite and keep a T169-like “fixed-magnitude + uniform g4 lift” shape,
does any strictly-defined reference magnitude lift (from a small fixed family of reference distributions) admit a sign solution?
```

So `P479`:

1. fixes a small finite family of reference distributions `q` on `Z_12`,
2. for each reference `q`, defines a magnitude lift `|vpsi_i| := sqrt(rho_*^2 * q_i)` from the strict QW-2122 scalars,
3. defines a uniform `g4` lift by `g4 := lambda_psi_strict / sum_i q_i^2`,
4. exhaustively scans all `2^11=2048` sign vectors (fixing `s0=+1`),
5. evaluates the `R18` declared `pair1` residual zero equations under the same conditional `N477` Yukawa-free diagonal residual rewrite.

## Inputs

1. QW-2122 strict scalar vacuum:
   - `report_qw2122_psi_potential_diagonal_floor_gate.json`
2. Strict-derived `alpha_geo`:
   - `fundamental_action_reconstruction/generated/alpha_geo_strict_derived_v1.json`
3. `R14` strict frozen kernel specialization (`K_total`):
   - `fundamental_action_reconstruction/generated/r14_explicit_frozen_kernel_specialization_packet_for_host_matching_route.json`
4. `R15` diagonal floor `m0^2`:
   - `fundamental_action_reconstruction/generated/r15_explicit_host_scalar_floor_embedding_packet_for_host_matching_route.json`
5. `R18` pair1 coefficient-class reduction and exact zero system:
   - `fundamental_action_reconstruction/generated/r18_pair1_residual_declared_pullback_coefficient_class_reduction_packet_for_host_matching_route.json`

## Output artifacts

- `fundamental_action_reconstruction/generated/pair1_residual_zero_equations_reference_magnitude_family_scan_under_n477_v1.json`
- `fundamental_action_reconstruction/generated/p479_current_strict_reference_magnitude_family_sign_scan_for_r18_pair1_zero_equations_probe.json`
- `fundamental_action_reconstruction/generated/p479_current_strict_reference_magnitude_family_sign_scan_for_r18_pair1_zero_equations_probe_summary.json`

## Hard limits (no false pass)

This probe does **not** claim:

1. any strict vacuum stationarity witness (it uses conditional `N477`),
2. any strict zero witness for the declared `pair1` residual equations (unless later promoted by a theorem-level discharge),
3. host matching (`C10_B1`) discharge,
4. selector-relevant physical canonicalization inside the `QW-2191` `O(2)` family,
5. strict-core selector closure / admissible `S_sel_int`,
6. global `QW-2191` discharge,
7. ToE closure.

