# P2676/S1626 Sigma->F301 pre-collapse naturality-square audit

Status: `P2676_PRECOLLAPSE_NATURALITY_SQUARE_AUDIT_NO_FALSE_PASS`

## Content-first repo grep
- `naturality_square_content`: `48` hits
- `sigma_to_f301_content`: `110` hits
- `precollapse_content`: `2720` hits
- `chart_label_retaining_content`: `333` hits
- `nonprojector_content`: `259` hits
- `nonconvention_content`: `8701` hits
- `selector_source_guard_content`: `15148` hits

## Finite Boolean naturality witness
Total maps checked: `16`.
Formal square candidates: `2`.
Formal candidate tables: `[6, 9]`.
Formal convention classes: `['xnor_reversal_orientation', 'xor_orientation']`.
Passing export gate count: `0`.

## Obstruction certificates
- `finite_formal_candidates_exist`: `True` — Boolean naturality-like maps exist as formal XOR/XNOR orientation choices.
- `orientation_reversal_pair_survives`: `True` — The formal candidates appear in a reversal pair, so finite naturality alone does not pick an internal orientation.
- `internal_orientation_source_absent`: `True` — No audited row has an internal selector/source bit that chooses XOR over XNOR without convention.
- `export_gate_zero`: `True` — No formal Boolean map is promoted to a current Sigma->F301 typed-arrow export.

## Verdict
P2676 builds the requested finite naturality-square witness for S6. Exhausting all 16 Boolean maps finds two formal chart-equivariant, source-sensitive, nonprojector candidates, but they are exactly an XOR/XNOR reversal pair. Since the repo does not export an internal orientation/source selecting one member of that pair, the passing export gate remains zero. Thus the naturality-square route does not currently supply the missing Sigma->F301 typed arrow.
Decision: `P2676_PRECOLLAPSE_NATURALITY_SQUARE_AUDIT__FORMAL_XOR_XNOR_PAIR_EXISTS_BUT_NO_INTERNAL_ORIENTATION_SOURCE`.

## Next honest step
The next honest move is not another same-shape Boolean/naturality lift. Either provide a new internal orientation source that breaks the XOR/XNOR reversal inside the Sigma->F301 square, or record a no-go for the current S6/O3 typed-seed route and return to a different bridge-completion/source class. O4/O5 and L_total reopening remain blocked.

## Negative exports
- `precollapse_naturality_square_exported`: `False`
- `sigma_to_f301_nonconvention_arrow_exported`: `False`
- `xor_xnor_orientation_selected_internally`: `False`
- `s6_sigma_to_f301_typed_arrow_exported`: `False`
- `o3_chart_sensitive_pair12_typed_seed_exported`: `False`
- `boundary_square_cycle_typed_arrow_exported`: `False`
- `sector_swap_sourced_invariant_exported`: `False`
- `q_w_2191_discharged`: `False`
- `role_bearing_ltotal_reenabled`: `False`
- `toe_closure_claimed`: `False`
