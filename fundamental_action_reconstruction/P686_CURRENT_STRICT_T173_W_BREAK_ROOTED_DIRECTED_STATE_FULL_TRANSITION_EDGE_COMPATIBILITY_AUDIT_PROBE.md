# P686 Current Strict `T173` `w_break`-Rooted Directed State Full Transition-Edge Compatibility Audit Probe (No False‑PASS)

Status: `P686_CURRENT_STRICT_T173_W_BREAK_ROOTED_DIRECTED_STATE_FULL_TRANSITION_EDGE_COMPATIBILITY_AUDIT_PROBE_SPEC_NO_FALSE_PASS`  
As of: `2026-03-17`

## Goal

After `F684/P684`, the repo exports one explicit **global directed (sign-sensitive) vector representative section**
on `C_v1` in a tracked convention scope:

- `SelectorState_global_C_v1_directed_rooted_transport_from_S_sel_int_w_break_strict_convention_v1`.

That object was constructed using **rooted** chart transports `O_1m` (from `pair1`) and therefore only directly audits
rooted alignment.

This probe audits the next honest compatibility question:

> Is the exported directed representative section compatible with the **full exported global transition/gluing object**
> on `C_v1` (`SelectorTransition_global_C_v1_strict_v1`) on every declared overlap edge `{pair1..pair5}`?

This remains **below** any strict physical sign datum:

- the global transition object is projector/axis-level (`alpha mod π` on several edges),
- `N512` forbids operator-level groupoid promotion,
- therefore any directed/vector-level consistency here is still a **section/convention audit** and does not imply
  kernel-alone/global `QW-2191` discharge.

## Inputs

1. Directed representative section:
   - `generated/selector_state_global_c_v1_directed_rooted_transport_from_s_sel_int_w_break_strict_convention_v1.json` (`F684`)
2. Global transition/gluing object on `C_v1`:
   - `generated/selector_transition_global_c_v1_strict_v1.json` (`F469`)
3. The per-edge operator files referenced by the global transition object (explicit `O_ij` matrices).

## Audit definition (edgewise)

For each declared transition edge `pairi_to_pairj` with operator `O_ij`:

1. compute `v := O_ij u_i` using the exported directed vectors `u_i` from `F684`,
2. test directed compatibility: `v ≈ u_j`,
3. also test the weaker projective compatibility: `v ≈ ± u_j` (line preservation).

The probe records:

- which edges preserve the line up to sign,
- which edges preserve the chosen directed sign,
- and where sign flips are present (still not a contradiction at projector level, but crucial for honest sign-datum discipline).

## Output

- `generated/p686_current_strict_t173_w_break_rooted_directed_state_full_transition_edge_compatibility_audit_probe.json`
- `generated/p686_current_strict_t173_w_break_rooted_directed_state_full_transition_edge_compatibility_audit_probe_summary.json`

## Hard limits (no false pass)

This probe does **not**:

- claim a directed/sign-sensitive physical orientation datum in strict core,
- claim kernel-alone/global `QW-2191` discharge,
- claim ToE closure,
- promote projector-level cocycle data into operator-level groupoid identities (`N512` boundary).

