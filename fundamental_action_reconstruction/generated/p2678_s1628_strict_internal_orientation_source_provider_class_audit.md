# P2678/S1628 strict internal orientation-source provider-class audit

Status: `P2678_STRICT_INTERNAL_ORIENTATION_SOURCE_PROVIDER_CLASS_AUDIT_NO_FALSE_PASS`

## Content-first repo grep
- `different_provider_class_content`: `1433` hits
- `oriented_torsor_content`: `34` hits
- `legacy_scalar_even_content`: `733` hits
- `symmetry_breaking_content`: `783` hits
- `sigma_f301_binding_content`: `105` hits
- `forbidden_collapse_content`: `10416` hits
- `closure_guard_content`: `15984` hits

## C2 equivariant provider enumeration
Sign maps checked: `4`.
Orientation-odd sign maps: `2`.
Formal odd provider classes: `['oriented_c2_torsor', 'spin_pin_orientation_source', 'boundary_symmetry_breaking_source']`.
Passing provider-class count: `0`.

## Provider rows
- `legacy_scalar_beta_tors` (even_scalar): formal_odd_domain=`False`, exported=`False`, binding=`False`, pass=`False`
- `axis_or_projector_quotient` (quotient_even): formal_odd_domain=`False`, exported=`False`, binding=`False`, pass=`False`
- `q_basis_terminal_branch` (postcollapse_choice): formal_odd_domain=`False`, exported=`False`, binding=`False`, pass=`False`
- `declared_orientation_convention` (fiat_choice): formal_odd_domain=`False`, exported=`False`, binding=`False`, pass=`False`
- `oriented_c2_torsor` (orientation_odd_torsor): formal_odd_domain=`True`, exported=`False`, binding=`False`, pass=`False`
- `spin_pin_orientation_source` (orientation_odd_spin_pin): formal_odd_domain=`True`, exported=`False`, binding=`False`, pass=`False`
- `boundary_symmetry_breaking_source` (orientation_odd_boundary): formal_odd_domain=`True`, exported=`False`, binding=`False`, pass=`False`

## Provider obligation lattice
Total states: `128`; passing states: `1`.
Current Hamming distance to pass: `5`.
Missing current obligations: `['orientation_odd_domain_exported', 'nonfiat_symmetry_breaking_or_torsor_law_exported', 'precollapse_sigma_side_binding_exported', 'surviving_f301_carrier_binding_exported', 'nonquotient_nonprojector_transport_exported']`.

## Verdict
P2678 moves to the different provider class requested after P2677.  The finite C2 enumeration shows exactly what kind of object could break XOR/XNOR without fiat: an orientation-odd torsor, spin/Pin orientation source, or boundary symmetry-breaking source.  Legacy scalar beta_tors, quotient/projector axes, Q_basis branches, and declared conventions cannot supply an odd strict sign.  The formal odd provider classes exist as admissible shapes, but none is currently exported together with a pre-collapse Sigma->F301 binding, so S6/O3 is not reopened.
Decision: `P2678_STRICT_INTERNAL_ORIENTATION_SOURCE_PROVIDER_CLASS_AUDIT__FORMAL_ODD_PROVIDERS_IDENTIFIED_BUT_NOT_EXPORTED`.

## Next honest step
The next honest proof-grade step is a construction-or-no-go audit for one concrete orientation-odd provider object, preferably the smallest C2 torsor/spin-boundary source with an explicit action on Sigma_sel_src_target_v1 and a typed arrow into the surviving F301 carrier.  If no such object can be exported, keep the S6/O3 route closed and shift the main bridge work back to the legacy->strict completion bridge source terms.

## Negative exports
- `strict_internal_orientation_source_theorem_exported`: `False`
- `oriented_torsor_provider_exported`: `False`
- `precollapse_sigma_f301_binding_exported`: `False`
- `xor_xnor_reversal_broken`: `False`
- `q_w_2191_discharged`: `False`
- `s6_reopened`: `False`
- `o3_reopened`: `False`
- `o4_o5_allowed`: `False`
- `role_bearing_ltotal_reenabled`: `False`
- `toe_closure_claimed`: `False`
