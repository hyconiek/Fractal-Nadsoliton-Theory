# P464 Current Strict Diagonal/Local `Psi` Hessian Eigenprojectors → Pair‑Plane Weight Profile Probe

Status: `P464_EXECUTED_EIGENPROJECTORS_TO_PAIR_PLANE_WEIGHT_PROFILE_PROBE_NO_FALSE_PASS`  
As of: `2026-03-15`

## Goal

`F459` exports a strict-derived **numeric** diagonal/local `Psi` Hessian instantiation `H_psi` together with an eigenbasis representative,
while `F460` exports the corresponding **sign‑gauge‑invariant** rank‑one spectral projectors:

```text
P_j := |v_j><v_j|,   j=0..11.
```

Independently, the diagonal/local lane exports a strict-derived mode-index assignment basis (`F453`) whose pair‑plane decomposition on `n=12` is:

```text
e0,  pair1, pair2, pair3, pair4, pair5,  e6.
```

This probe exports a **pair‑plane weight profile** for each `H_psi` eigenprojector:

```text
w_{j,label} := tr( Π_label · P_j )
```

where `Π_label` is the orthogonal projector onto the corresponding pair plane (or rank‑one subspace for `e0`/`e6`).

This construction is:
- sign‑gauge‑invariant (uses projectors),
- axis‑invariant inside each pair plane (uses the full plane projector, not a chosen axis),
- lane‑scoped (depends on the exported strict-derived diagonal/local objects).

It does **not** claim any global selector transition/gluing object and does **not** claim ToE closure.

## Inputs

1. `F460`
   - `fundamental_action_reconstruction/generated/psi_hessian_diagonal_local_strict_derived_eigenprojectors_v1.json`
2. `F453`
   - `fundamental_action_reconstruction/generated/mode_index_assignment_canonical_local_diagonal_strict_derived_v1.json`

## Output artifacts

- `fundamental_action_reconstruction/generated/p464_current_strict_diagonal_local_psi_hessian_eigenprojectors_to_pair_plane_weight_profile_probe.json`
- `fundamental_action_reconstruction/generated/p464_current_strict_diagonal_local_psi_hessian_eigenprojectors_to_pair_plane_weight_profile_probe_summary.json`

## Hard limits (no false pass)

This probe does **not** claim:

1. that the `F453` mode-index assignment basis diagonalizes the full `H_psi`,
2. any global selector transition/gluing object (`H40_B1` remains),
3. any discharge of `QW-2191` beyond the declared lane-scoped instantiations,
4. strict-core selector closure / admissible `S_sel_int`,
5. ToE closure.

