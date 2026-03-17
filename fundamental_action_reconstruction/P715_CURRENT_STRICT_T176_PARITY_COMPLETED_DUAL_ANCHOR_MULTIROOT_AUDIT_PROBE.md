# P715 Current Strict `T176` Parity-Completed Dual-Anchor Multi-Root Audit Probe (No False-PASS)

Status: `P715_CURRENT_STRICT_T176_PARITY_COMPLETED_DUAL_ANCHOR_MULTIROOT_AUDIT_PROBE_SPEC_NO_FALSE_PASS`
As of: `2026-03-17`

## Goal

After `P714`, the next honest question is no longer merely:

```text
do we need an even-or-mixed parity anchor component?
```

but the sharper continuation:

```text
what happens if we use the minimal already-exported parity-completed pair
(w_break, w_ref_unnormalized)
under one fixed slot-free root rule?
```

This probe audits that exact continuation without introducing any new provider
object or hidden selector slot.

## Inputs

1. `generated/f647_strict_witness_provider_export_packet_for_seed_v1_realization_attempt.json`
   - exported strict-core payloads `w_break_by_x` and `w_ref_unnormalized_by_x`.
2. `generated/a_m_pairm_orientation_projector_operator_strict_core_v1.json` for `m=1..5`
   - exported representative vectors `u_m` and chart coordinates.
3. `generated/selector_atlas_pair12345_oriented_transport_mod_2pi_strict_convention_v1.json`
   - exported oriented transport family on `{pair1..pair5}`.
4. `generated/selector_output_operator_global_c_v1_seed_v1_promoted_strict_v1.json`
   - fixed output-basis readout used only to compare rooted sections.
5. `generated/p714_current_strict_t176_w_break_parity_root_support_profile_audit_probe_summary.json`
   - current parity-corridor boundary for the odd-only anchor class.

## Audit rule

For each possible root chart `pair_r`, define the strict-core root anchor rule:

```text
s_r := sign( Σ_x w_break(x) u_r(x) )          if that scalar is nonzero,
     := sign( Σ_x w_ref(x)   u_r(x) )          otherwise, if that scalar is nonzero,
     := undefined                              otherwise.
```

where `w_ref := w_ref_unnormalized_by_x` is used exactly as exported by `F647`.

Then:

1. transport `s_r u_r` through the exported oriented operator family,
2. align with each target representative `u_m`,
3. compare the recovered chartwise sign vectors and output vectors across all
   roots:
   - first **exactly**,
   - then only **up to one global sign**.

## Result discipline

The strongest positive result would be:

1. all roots are supported,
2. all roots recover the same exact directed section,
3. and the result still remains below `T176` because the transport family is
   convention-layer only.

If exact agreement fails, the probe must still report the strongest surviving
subcase honestly. In particular, if all roots are supported but agree only up
to one global `Z2` sign, that is a **projective/root-orbit** survival result,
not a directed-sign discharge.

## Hard limits

This probe does **not**:

- claim a strict-core directed/sign-sensitive physical orientation datum,
- claim kernel-alone/global `QW-2191` discharge,
- claim `T176` discharge,
- promote convention-layer oriented transport into a strict-core/global
  provider,
- claim ToE closure.
