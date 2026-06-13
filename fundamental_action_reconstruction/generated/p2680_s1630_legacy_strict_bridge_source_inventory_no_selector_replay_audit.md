# P2680/S1630 legacy->strict bridge-source inventory without selector replay

Status: `P2680_LEGACY_STRICT_BRIDGE_SOURCE_INVENTORY_NO_SELECTOR_REPLAY_AUDIT_NO_FALSE_PASS`

## Content-first repo grep
- `legacy_strict_kernel_content`: `2533` hits
- `amplitude_normalization_content`: `8037` hits
- `damping_compression_content`: `1012` hits
- `positive_beta_source_content`: `460` hits
- `canonical_unit_content`: `386` hits
- `forbidden_replay_content`: `16391` hits

## Source atom inventory
- `legacy_kernel_source_layer_visible` (input discipline): formal=`True`, source_export=`True`, role_safe=`True` — K_legacy_ont and canonical D_f/alpha_geo/beta_tors layer remain visible as bridge input, not final strict identity.
- `strict_kernel_target_layer_visible` (target discipline): formal=`True`, source_export=`True`, role_safe=`True` — K_strict_gate is available as operational strict target, not silently identical to the legacy kernel.
- `alpha_geo_scalar_shape_normalization` (amplitude/normalization): formal=`True`, source_export=`True`, role_safe=`False` — Scalar shape normalization exists, but role-safe amplitude absorption/physical-role transfer is not exported.
- `amplitude_role_safe_source` (amplitude/normalization): formal=`False`, source_export=`False`, role_safe=`False` — No theorem says legacy alpha_geo roles survive or transform into strict Lagrangian/value roles.
- `fractal_pushforward_linear_to_power_damping_shape` (damping/compression): formal=`True`, source_export=`True`, role_safe=`False` — The formal q(d)=d^(9/5) shape route can change linear damping form into power damping shape, but it is not the full coefficient source.
- `target_independent_positive_beta_or_z_beta_source` (damping/compression): formal=`False`, source_export=`False`, role_safe=`False` — The missing source atom is a target-independent positive beta/Z_beta coefficient source, not a normalization fit.
- `canonical_length_or_uv_unit_source` (damping/compression): formal=`False`, source_export=`False`, role_safe=`False` — The beta=1 normalization orbit is understood, but no canonical length/UV unit theorem fixes the gauge physically.
- `selector_phase_orientation_source` (phase/topological selector): formal=`False`, source_export=`False`, role_safe=`False` — Deliberately not reopened by P2680; P2679 forbids selector replay without a new object.

## Bridge-source lattice
Total states: `256`; passing states: `1`.
Current Hamming distance to bridge-source pass: `3`.
Missing current obligations: `['amplitude_role_safe_source_exported', 'positive_beta_z_beta_source_exported', 'canonical_length_uv_unit_exported']`.

## Component readiness
- `amplitude/normalization`: ready_now=`False`; missing=['alpha_geo_scalar_shape_normalization', 'amplitude_role_safe_source'] — one or more source atoms remain missing or not role-safe
- `damping/compression`: ready_now=`False`; missing=['fractal_pushforward_linear_to_power_damping_shape', 'target_independent_positive_beta_or_z_beta_source', 'canonical_length_or_uv_unit_source'] — one or more source atoms remain missing or not role-safe
- `phase/topological selector`: ready_now=`False`; missing=['selector_phase_orientation_source'] — one or more source atoms remain missing or not role-safe
- `role transfer`: ready_now=`False`; missing=['full bridge completion', 'claim-by-claim role-transfer theorem'] — Role transfer is downstream of full bridge completion and remains forbidden by guardrail.

## Verdict
P2680 follows the P2679 pivot and does not reopen selector, tau_src->pair12, or beta_tors->chi11.  The bridge-source inventory finds real formal material: legacy/strict kernel layers are visible, scalar alpha_geo shape normalization exists, and a fractal pushforward can supply the damping power-shape route.  The proof-grade bridge still fails because three non-selector source atoms are missing or not role-safe: amplitude role-safe source, target-independent positive beta/Z_beta source, and canonical length/UV unit source.
Decision: `P2680_LEGACY_STRICT_BRIDGE_SOURCE_INVENTORY__AMPLITUDE_AND_DAMPING_SOURCES_STILL_INCOMPLETE_NO_SELECTOR_REPLAY`.

## Next honest step
The next honest proof-grade move is not selector work.  Choose one missing non-selector atom and run a construction-or-no-go audit, with highest leverage on the target-independent positive beta/Z_beta source or the canonical length/UV unit source.  Only after those source atoms are exported should role-transfer auditing be attempted.

## Negative exports
- `full_legacy_strict_bridge_completed`: `False`
- `amplitude_role_safe_source_exported`: `False`
- `positive_beta_renormalization_source_exported`: `False`
- `canonical_length_uv_unit_exported`: `False`
- `selector_replay_used`: `False`
- `beta_tors_chi11_reopened`: `False`
- `role_transfer_theorem_exported`: `False`
- `role_bearing_ltotal_reenabled`: `False`
- `q_w_2191_discharged`: `False`
- `toe_closure_claimed`: `False`
