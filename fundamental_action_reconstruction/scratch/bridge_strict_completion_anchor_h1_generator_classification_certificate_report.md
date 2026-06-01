# Anchor vs H1 generator classification certificate

Status: `anchor-is-c0-gauge-fix-not-h1-generator__closed-cycle-generator-is-odd-parity-c1-class`

## Result

The d=0 left anchor is classified as a C0 gauge-fixing datum, not as the
closed-cycle H1 generator.  The closed-cycle generator is represented by an
odd-parity C1 edge cochain such as the single closing edge `11->0`.

## Summary

- `path_h1_zero`: `True`
- `cycle_h1_one`: `True`
- `cycle_image_equals_even_parity_kernel`: `True`
- `closing_edge_odd_parity_generator`: `True`
- `alternate_odd_edge_same_generator_class`: `True`
- `audited_closed_cochain_is_exact_zero_h1_class`: `True`
- `left_anchor_is_c0_gauge_fix_not_c1_generator`: `True`
- `ai_opinion_as_stated_rejected`: `True`
- `ai_opinion_weak_form_accepted`: `True`
- `selector_source_remains_open`: `True`
- `gf2_path_solution_unique_inherited`: `True`

## Opinion audit

- `accepted_for_closed_cycle_only`: H^1(Z12;GF(2)) is one-dimensional on the closed 12-cycle. (cycle report has rank(delta)=11, dim C1=12, so dim H^1=1)
- `rejected_type_mismatch`: The d=0 left anchor is itself the H^1 generator. (b(0) is a C0 gauge-fixing value; an H^1 generator is a C1/im(delta) edge-cochain class with odd cycle parity)
- `rejected_for_audited_zero_closure`: The audited phase pattern represents the nonzero H^1 generator. (audited path plus forced closing bit has total parity 0 and is exact on the closed cycle)
- `accepted`: The missing selector/source is still open. (the generator classification gives the address of the obstruction, not a source theorem or QW-2191 discharge)

## Hard limits

- No beta_tors -> chi_11 theorem is claimed.
- No strict-core selector source is exported.
- No QW-2191 selector discharge is claimed.
- No legacy physical-role transfer to K_strict_gate is claimed.
- No unqualified identity K_legacy_ont == K_strict_gate is claimed.
- No ToE closure is claimed.
