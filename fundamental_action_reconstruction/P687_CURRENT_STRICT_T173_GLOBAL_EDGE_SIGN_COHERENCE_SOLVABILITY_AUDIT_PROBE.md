# P687 Current Strict `T173` Global Edge Sign-Coherence Solvability Audit Probe (No False‑PASS)

Status: `P687_CURRENT_STRICT_T173_GLOBAL_EDGE_SIGN_COHERENCE_SOLVABILITY_AUDIT_PROBE_SPEC_NO_FALSE_PASS`  
As of: `2026-03-17`

## Goal

`P686` audits the exported `w_break`-rooted directed selector state representative (`F684`) against the **full** exported global transition/gluing object
on `C_v1` (`F469`). It confirms edgewise compatibility only **up to sign** and records which overlap edges force a sign flip under the currently exported
axis-only (`α mod π`) transport representatives.

This probe sharpens the next honest question:

> Can the recorded edgewise sign flips be removed purely by a **per-chart sign relift** (`u_i -> t_i u_i`, `t_i∈{±1}`),
> while keeping the exported transition operators fixed?
>
> Equivalently: does there exist a `Z2`-valued 0‑cochain `{t_i}` such that for every exported overlap edge `i→j`,
> the transported directed vectors satisfy `O_ij (t_i u_i) ≈ (t_j u_j)` (no sign flip on any edge)?

If the answer is **no**, this identifies a genuine **global edge sign-coherence obstruction** relative to the fixed exported axis-only transition representatives:
the sign-flip pattern is not eliminable by chart sign choices alone.

This remains strictly below any physical sign datum claim; it is a compatibility/obstruction statement about exported artifacts.

## Inputs

1. `P686` full-edge directed compatibility audit artifact:
   - `generated/p686_current_strict_t173_w_break_rooted_directed_state_full_transition_edge_compatibility_audit_probe.json`

## Output

- `generated/p687_current_strict_t173_global_edge_sign_coherence_solvability_audit_probe.json`
- `generated/p687_current_strict_t173_global_edge_sign_coherence_solvability_audit_probe_summary.json`

## Hard limits (no false pass)

This probe does **not**:

- claim a directed/sign-sensitive physical orientation datum in strict core,
- claim kernel-alone/global `QW-2191` discharge,
- claim ToE closure.

