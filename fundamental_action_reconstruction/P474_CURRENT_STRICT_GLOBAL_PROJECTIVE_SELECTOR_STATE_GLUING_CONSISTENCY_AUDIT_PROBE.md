# P474 Current Strict Global Projective Selector State — Gluing/Transport Consistency Audit Probe (No False‑PASS)

Status: `P474_EXECUTED_CURRENT_STRICT_GLOBAL_PROJECTIVE_SELECTOR_STATE_GLUING_CONSISTENCY_AUDIT_PROBE_NO_FALSE_PASS`  
As of: `2026-03-16`

## Goal

After `F469/N515` and `F470/N516`, strict core exports on the declared strict domain `C_v1`:

1. a **global selector atlas object** `SelectorAtlas_global_C_v1_strict_v1`,
2. a **global selector transition/gluing object** `SelectorTransition_global_C_v1_strict_v1`,
3. a **global projective selector state object** `SelectorState_global_C_v1_projective_strict_v1`
   (ray/projector semantics; residual sign is gauge at state level).

This probe audits one no‑false‑pass hygiene condition:

```text
the exported global projective selector state is actually glued/transported
consistently by the exported global selector transition operators, at
projector (ray) level.
```

## What it checks (declared scope)

1. Load the global projective selector state and its chart-local projectors `A_m(pair_m)` (`m=1..5`).
2. Load the global selector transition operators `O_{ij}` exported by `SelectorTransition_global_C_v1_strict_v1`.
3. For each directed edge `i→j`, check:
   - orthogonality: `O_{ij}^T O_{ij} ≈ I`,
   - projector transport: `P_j ≈ O_{ij} P_i O_{ij}^T`, where `P_m := |u_m><u_m|` is reconstructed from the exported `u_m`.
4. For each ordered triple `(i,j,k)` where edges exist, check projector‑level cocycle on the exported projector section:
   `O_{jk} O_{ij} P_i (O_{jk} O_{ij})^T ≈ O_{ik} P_i O_{ik}^T`.

This is intentionally **projector‑level** (ray‑level) only.
It does not attempt to upgrade residual sign into a strict sign‑sensitive physical datum.

## Exported artifacts

Executed by:

```text
python3 fundamental_action_reconstruction/p474_current_strict_global_projective_selector_state_gluing_consistency_audit_probe.py
```

Exports:

- `fundamental_action_reconstruction/generated/p474_current_strict_global_projective_selector_state_gluing_consistency_audit_probe.json`
- `fundamental_action_reconstruction/generated/p474_current_strict_global_projective_selector_state_gluing_consistency_audit_probe_summary.json`

## Hard limits (no false pass)

`P474` does **not** claim:

1. a sign-sensitive / directed selector state datum (lifting residual `Z2`),
2. strict-core selector closure / admissible `S_sel_int`,
3. global discharge of `QW-2191`,
4. ToE closure.

