# P683 Current Strict `T173` Rooted-Transport Sign-Lift From `S_sel_int` `w_break` Audit Probe

Status: `P683_CURRENT_STRICT_T173_ROOTED_TRANSPORT_SIGN_LIFT_FROM_S_SEL_INT_W_BREAK_AUDIT_PROBE_SPEC_NO_FALSE_PASS`  
As of: `2026-03-17`

## Goal

`P681` audited a strict-core-only attempt to lift residual `Z2` sign on all charts using only the exported seed payload weights from `F647` and per-chart dot products.
That candidate yielded an output-sign mismatch and was packaged as a strict boundary (`N681`).

This probe audits a different, still minimal, **rooted** candidate:

> use the exported reflection-breaking seed weight `w_break` **only to fix sign on the seed/root chart** (`pair1`),
> then propagate sign to other charts via the exported chart-transport operators `O_1m` by requiring vector-level transport
> `O_1m u_1^* ≈ u_m^*` (choosing the sign of the representative `u_m` accordingly).

This is intended as a narrow “does the currently exported transition infrastructure support a consistent *section choice* once a root sign is fixed?” audit,
without claiming any strict physical sign datum.

## Inputs

1. `generated/f647_strict_witness_provider_export_packet_for_seed_v1_realization_attempt.json`
   - provides the exported strict seed payload reflection-breaking weights `w_break_by_x`.
2. `generated/a_m_pairm_orientation_projector_operator_strict_core_v1.json` for `m=1..5`
   - provides exported representative vectors `u_m` and their chart coordinates.
3. Transport operators rooted at `pair1`:
   - `generated/o12_pair1_pair2_selector_chart_transport_operator_strict_derived_from_sigma_int_alpha12_v1.json`
   - `generated/o13_pair1_pair3_selector_chart_transport_operator_axis_only_alpha13_mod_pi_strict_core_v1.json`
   - `generated/o14_pair1_pair4_selector_chart_transport_operator_axis_only_alpha14_mod_pi_strict_core_v1.json`
   - `generated/o15_pair1_pair5_selector_chart_transport_operator_axis_only_alpha15_mod_pi_strict_core_v1.json`
4. `generated/selector_output_operator_global_c_v1_seed_v1_promoted_strict_v1.json`
   - chartwise output channels `Y_sel(pair_m)` used only to audit output-sign consistency after applying the rooted sign lift.

## Rooted sign-lift rule (audit definition)

1. Root sign on `pair1`:

   ```text
   s_1 := sign( Σ_x w_break(x) u_1(x) )   (if nonzero)
   u_1^* := s_1 u_1
   ```

2. Propagate to each `pair_m` (m=2..5) using the exported `O_1m` operator:

   ```text
   v_m := O_1m u_1^*
   s_m := sign( <v_m, u_m> )
   u_m^* := s_m u_m
   ```

3. (sanity) apply the sign to chart coordinates and push through `Y_sel(pair_m)`; record the sign of the `o_+` component.

## Acceptance (no false pass)

The probe reports:

- whether this rooted rule supplies a sign lift on all charts, and
- whether the resulting `o_+` output sign is consistent across `{pair1..pair5}`.

It must keep explicit that this is still not a strict physical orientation datum:

- the rooted propagation uses axis-only transport edges on some overlaps (`alpha mod π`), which are only projectively meaningful,
- therefore any oriented-vector-level lift is a *section/convention choice* unless separately upgraded to an oriented `mod 2π` transport class.

## Hard limits

This probe does **not**:

- claim a directed/sign-sensitive physical orientation datum in strict core,
- claim kernel-alone/global `QW-2191` discharge,
- claim ToE closure,
- promote projector-level cocycle into an operator-level groupoid (`N512` boundary).

