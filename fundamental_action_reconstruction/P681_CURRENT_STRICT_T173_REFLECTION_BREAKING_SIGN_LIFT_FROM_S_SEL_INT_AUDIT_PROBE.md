# P681 Current Strict T173 Reflection-Breaking Sign-Lift From `S_sel_int` Audit Probe

Status: `P681_CURRENT_STRICT_T173_REFLECTION_BREAKING_SIGN_LIFT_FROM_S_SEL_INT_AUDIT_PROBE_SPEC_NO_FALSE_PASS`  
As of: `2026-03-17`

## Goal

After projective strict-core selector closure is discharged (`N680`), the remaining strict frontier tracked by `T173/N679` includes:

1. kernel-alone/global `QW-2191` discharge (still unclaimed), and
2. any **directed/sign-sensitive physical orientation datum** in strict core beyond projective/ray semantics.

This probe attacks only (2), in the narrowest way:

> Audit whether the currently exported strict-core internal selector-source seed object witness (`S_sel_int_strict_core_source_object_v1` export witness provider `F647`)
> already contains explicit site-indexed weight ingredients that can deterministically lift the residual `Z2` sign of the
> exported representative vectors on `{pair1..pair5}` in a single global rule:
>
> - even reference weights (`w_ref`) can, in principle, distinguish sign on **even** axes (cos-like representatives),
> - reflection-breaking weights (`w_break`) are required to distinguish sign on **odd** axes (sin-like representatives).

In other words: can the exported strict-core reflection-breaking ingredient distinguish `u` from `-u` beyond `pair1` without adding any additional premise?

## Inputs

1. `generated/f647_strict_witness_provider_export_packet_for_seed_v1_realization_attempt.json`
   - contains the strict-core source-object candidate payload including a reflection-breaking weight array `w_break_by_x`.
2. `generated/a_m_pairm_orientation_projector_operator_strict_core_v1.json` for `m=1..5`
   - contains the exported representative vectors `u_m` (projector-level objects but with explicit representative vectors recorded) and their chart coordinates.
3. `generated/selector_output_operator_global_c_v1_seed_v1_promoted_strict_v1.json`
   - provides the exported chartwise output channels `Y_sel(pair_m)` (used only for a sanity check of output-sign consistency after any proposed sign lift).

## Acceptance (no false pass)

For each `pair_m (m=1..5)`:

1. compute the scalars:
   - `dot_ref_m := Σ_x w_ref(x) u_m(x)`,
   - `dot_break_m := Σ_x w_break(x) u_m(x)`,
2. define a deterministic sign-lift candidate by the first available nonzero scalar:
   - if `|dot_ref_m| > tol`, take `s_m := sign(dot_ref_m)`,
   - else if `|dot_break_m| > tol`, take `s_m := sign(dot_break_m)`,
   - else: no sign-lift candidate on that chart from these exported weights,
3. define the oriented representative `u_m^* := s_m u_m` where available,
4. (sanity) push the oriented chart coordinates through the exported output channel `Y_sel(pair_m)` and record the sign of the `o_+` component.

The probe must report:

- whether a single strict-core reflection-breaking weight ingredient suffices to lift sign on **all** charts,
- or only on a strict subset (partial),
- or on none (fail),

without claiming that any such lift is physically canonical (that would require additional theorem-level work).

## Hard limits

This probe does **not**:

- discharge any directed/sign-sensitive physical orientation datum in strict core,
- claim kernel-alone/global `QW-2191` discharge,
- claim ToE closure,
- upgrade any premise-based directed objects into strict-core directed objects.
