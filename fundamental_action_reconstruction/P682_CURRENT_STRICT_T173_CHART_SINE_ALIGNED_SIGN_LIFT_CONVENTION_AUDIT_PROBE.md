# P682 Current Strict `T173` Chart-Sine-Aligned Sign-Lift Convention Audit Probe

Status: `P682_CURRENT_STRICT_T173_CHART_SINE_ALIGNED_SIGN_LIFT_CONVENTION_AUDIT_PROBE_SPEC_NO_FALSE_PASS`  
As of: `2026-03-17`

## Goal

After projective strict-core selector closure is discharged (`N680`), the remaining strict frontier tracked by `T173/N679` still includes any **directed/sign-sensitive**
global closure outcome in strict core beyond projective/ray semantics.

`P681` audited a narrow candidate:

> use only the exported strict-core seed payload weights (`F647`) to deterministically lift residual `Z2` sign chartwise.

It reported a directed output sign mismatch under that strict-core-only candidate.

This probe performs a different, explicitly **non-physical convention-layer** audit:

> if one allows a *chart-sine-aligned* reflection-breaking weight family
> `w_break_m(x) := sigma_int * w_ref(x) * s_m(x)` (where `s_m` is the declared chart sine basis on `Z_12`),
> does this supply a deterministic per-chart sign lift that yields a globally consistent output sign for the exported global output channels?

This probe is intended to separate two statements cleanly:

1. strict-core-only sign lift (still obstructed by `P681` / packaged by `N681`), vs
2. a deterministic **chart-convention** sign lift that can be used where sign is only a tracked gauge (not a physical datum).

## Inputs

1. `generated/f647_strict_witness_provider_export_packet_for_seed_v1_realization_attempt.json`
   - provides the exported strict seed payload **even** reference weights `w_ref`.
2. `generated/theta_pair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1.json`
   - provides the strict-derived `sigma_int` value used in the sigma-int corridor.
3. `generated/a_m_pairm_orientation_projector_operator_strict_core_v1.json` for `m=1..5`
   - provides the exported representative vectors `u_m` and their chart coordinates.
4. `generated/selector_output_operator_global_c_v1_seed_v1_promoted_strict_v1.json`
   - provides exported chartwise output channels `Y_sel(pair_m)` to audit output-sign consistency after applying the sign-lift convention.

## Convention-level sign-lift rule (audit definition)

For each `pair_m (m=1..5)`:

1. compute `dot_ref_m := Σ_x w_ref(x) u_m(x)`,
2. define the chart-sine-aligned weight:

   ```text
   w_break_m(x) := sigma_int * w_ref(x) * s_m(x),
   s_m(x) := sin(2π m x / 12)
   ```

3. compute `dot_break_m := Σ_x w_break_m(x) u_m(x)`,
4. propose a deterministic sign lift:
   - if `|dot_ref_m| > tol`, take `s_m := sign(dot_ref_m)`,
   - else if `|dot_break_m| > tol`, take `s_m := sign(dot_break_m)`,
   - else: no sign lift on that chart under this convention rule,
5. apply the sign to the exported chart coordinates and push through `Y_sel(pair_m)` as a sanity check:
   - report the sign of the `o_+` component.

## Acceptance (no false pass)

This probe reports only whether the above **convention rule** yields:

- a sign lift on all charts, and
- a globally consistent `o_+` output sign across charts.

It must also keep explicit that this is **not** a strict physical sign datum:

- the use of `s_m(x) = sin(2π m x/12)` depends on a non-`Aut(Z_12)`-invariant chart embedding convention of `Z_12` into `U(1)`,
- therefore any success here counts only as a *tracked convention layer*, not as a strict-core discharge of directed physical orientation.

## Hard limits

This probe does **not**:

- discharge any directed/sign-sensitive physical orientation datum in strict core,
- claim kernel-alone/global `QW-2191` discharge,
- claim ToE closure,
- upgrade any premise-based directed objects into strict-core directed objects.

