# P2620/S1570 P2618/P2619 bridge two-obstruction cut

Status: `P2620_TWO_OBSTRUCTION_BRIDGE_SOURCE_CUT_EXACT_TRUTH_TABLE_NO_FULL_BRIDGE_NO_LTOTAL_NO_QW2191_NO_TOE`

## Theorem

After P2618/P2619, any non-role-bearing legacy-to-strict bridge-source repair must supply both an independent nonlinear damping completion source and an independent orientation-odd selector source. Either one alone leaves a named obstruction alive.

## Proof

- P2618 proves that the strict exponent source eta=9/5 does not provide an exact scalar beta_tors -> beta damping completion map.
- P2619 proves that orientation-invariant legacy scalar or axis-only data cannot provide a C2-equivariant strict odd phase-sign selector.
- The damping obstruction and selector obstruction live in different output coordinates of the completion map: denominator compression versus phase/topological sign.
- Therefore repairing only the damping coordinate leaves the selector obstruction, and repairing only the selector coordinate leaves the damping obstruction.
- The bridge-source cut is the conjunction of the two independent repair atoms; the finite truth table has exactly one accepting row.

## Computed truth table

- `{'nonlinear_damping_completion_source': False, 'orientation_odd_selector_source': False}` -> bridge-source cut repaired `False`, missing `['nonlinear_damping_completion_source', 'orientation_odd_selector_source']`.
- `{'nonlinear_damping_completion_source': False, 'orientation_odd_selector_source': True}` -> bridge-source cut repaired `False`, missing `['nonlinear_damping_completion_source']`.
- `{'nonlinear_damping_completion_source': True, 'orientation_odd_selector_source': False}` -> bridge-source cut repaired `False`, missing `['orientation_odd_selector_source']`.
- `{'nonlinear_damping_completion_source': True, 'orientation_odd_selector_source': True}` -> bridge-source cut repaired `True`, missing `[]`.

## Shortcut rejection matrix

- `eta_9_5_exponent_source_only`: passes `False` — P2618 retains the exponent source but blocks exact scalar damping completion and selector source.
- `beta_tors_scalar_renormalization`: passes `False` — P2618 derivative obstruction blocks scalar damping completion; P2619 treats beta_tors as orientation-invariant.
- `axis_only_selector_up_to_Z2`: passes `False` — Axis-only data can reduce continuous degeneracy but leaves the strict odd sign unresolved.
- `GF2_bookkeeping_rank_or_cycle_data`: passes `False` — Combinatorial bookkeeping may classify constraints but is not a physical orientation/source premise under P2612/P2619.
- `damping_source_plus_no_selector`: passes `False` — Even a future damping completion source alone leaves the P2619 C2 selector obstruction.
- `selector_source_plus_no_damping_completion`: passes `False` — Even a future selector source alone leaves the P2618 nonlinear damping completion obstruction.
- `both_independent_sources_supplied`: passes `True` — This is the minimal bridge-source cut repair, still not a role-transfer theorem.

## Recommended next honest step

Pick one real source target, not a shortcut: either derive a nonlinear damping completion source beyond scalar beta_tors renormalization, or derive an orientation-odd selector source; full bridge-source repair needs both before role transfer can be rerun.

## Scope guards

No full bridge revalidation, no role-transfer revalidation, no role-bearing `L_total`, no `beta_tors -> chi11` reopening, no GF(2) bridge revalidation, no `QW-2191` discharge, and no ToE closure are exported.

## Fingerprint

`bd73b6a915764cb118b518dd8699d9b6b58468ce1f5a5852e1127ecbe0299860`
