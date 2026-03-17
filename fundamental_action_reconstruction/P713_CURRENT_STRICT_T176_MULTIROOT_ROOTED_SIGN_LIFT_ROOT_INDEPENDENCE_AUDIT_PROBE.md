# P713 Current Strict `T176` Multi-Root Rooted Sign-Lift Root-Independence Audit Probe (No False-PASS)

Status: `P713_CURRENT_STRICT_T176_MULTIROOT_ROOTED_SIGN_LIFT_ROOT_INDEPENDENCE_AUDIT_PROBE_SPEC_NO_FALSE_PASS`
As of: `2026-03-17`

## Goal

After `P712/N708`, the next honest continuation is not to pretend that a strict-core global provider already exists, but to ask a sharper question:

```text
is the currently surviving rooted directed-section route merely a pair1-root artifact,
or does it recover the same directed section from every chart root once the root sign is fixed locally by w_break
and propagated through the exported oriented transport family?
```

This probe audits that question on the current exported artifacts.

## Inputs

1. `generated/f647_strict_witness_provider_export_packet_for_seed_v1_realization_attempt.json`
   - exported strict-core reflection-breaking payload `w_break_by_x`.
2. `generated/a_m_pairm_orientation_projector_operator_strict_core_v1.json` for `m=1..5`
   - exported representative vectors `u_m` and chart coordinates.
3. `generated/selector_atlas_pair12345_oriented_transport_mod_2pi_strict_convention_v1.json`
   - exported lane-scoped oriented transport family on `{pair1..pair5}`.
4. `generated/selector_output_operator_global_c_v1_seed_v1_promoted_strict_v1.json`
   - fixed output-basis readout used only to compare rooted sections.
5. Existing rooted convention-layer result:
   - `generated/f684_first_exported_selector_state_global_c_v1_directed_rooted_transport_from_s_sel_int_w_break_packet_summary.json`

## Audit rule

For each possible root chart `pair_r`, define:

1. root sign from exported strict-core payload:

   ```text
   s_r := sign( Σ_x w_break(x) u_r(x) )
   ```

2. for every other chart `pair_m`, propagate sign by transporting `s_r u_r`
   through the exported oriented operator on the edge `{pair_r,pair_m}` and aligning with `u_m`.

The probe then tests whether the resulting sign vectors and output vectors are identical for all supported roots.

## Result discipline

The strongest positive result would be:

1. every root chart supports a nonzero `w_break` sign anchor,
2. transport alignment succeeds for every root-to-chart propagation,
3. all roots recover the same chartwise sign vector and output-space vector section,
4. the result still stays below `T176` because the transport family remains convention-layer only.

If that stronger result fails, the probe must still report the strongest surviving subcase honestly — for example a supported-root corridor on which the same convention-layer section is recovered, without promoting that to a global provider.

## Hard limits

This probe does **not**:

- claim a strict-core directed/sign-sensitive physical orientation datum,
- claim kernel-alone/global `QW-2191` discharge,
- claim that `T176` is discharged,
- promote convention-layer oriented transport into a strict-core/global provider,
- claim ToE closure.
