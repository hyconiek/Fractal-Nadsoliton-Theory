# Affine subgroup lattice chi_11 source obstruction

Status: `only-cyclic-and-D12-quotients-admit-chi11-unit5-or-unit7-kill-it`

## Finite model

- Ring: `Z_12`
- Enumerated supports: `792`
- Affine unit subgroup count: `5`
- All automorphism units: `[1, 5, 7, 11]`

## Lattice summary

- Admitting subgroup count: `2`
- Killing subgroup count: `3`
- Admitting subgroups: `['T_semidirect_{1}_cyclic_translation', 'T_semidirect_{1,11}_D12']`
- Killing subgroups: `['T_semidirect_{1,5}_unit5_axis', 'T_semidirect_{1,7}_unit7_axis', 'T_semidirect_{1,5,7,11}_full_affine_Aut']`
- Minimal nontrivial admitting reflection subgroup: `T_semidirect_{1,11}_D12`
- Minimal killing unit additions: `['unit5', 'unit7']`
- Full Aut admits chi_11: `False`

## Subgroup rows

- `T_semidirect_{1}_cyclic_translation` units `[1]`: orbits `66`, branch chi_11 well-defined `True`, classification `admits_branch_chi11_quotient`
- `T_semidirect_{1,11}_D12` units `[1, 11]`: orbits `38`, branch chi_11 well-defined `True`, classification `admits_branch_chi11_quotient`
- `T_semidirect_{1,5}_unit5_axis` units `[1, 5]`: orbits `38`, branch chi_11 well-defined `False`, classification `kills_branch_chi11_by_mixing_opposite_signs`
- `T_semidirect_{1,7}_unit7_axis` units `[1, 7]`: orbits `40`, branch chi_11 well-defined `False`, classification `kills_branch_chi11_by_mixing_opposite_signs`
- `T_semidirect_{1,5,7,11}_full_affine_Aut` units `[1, 5, 7, 11]`: orbits `25`, branch chi_11 well-defined `False`, classification `kills_branch_chi11_by_mixing_opposite_signs`

## Proof certificate

- `finite_domain`: All C(12,5)=792 supports are enumerated for every affine unit subgroup between translations and full Aut.
- `well_definedness_test`: A quotient can host branch chi_11 exactly when no quotient orbit contains both signs A1/A11=-1 and A5/A7=+1.
- `admitting_subgroups`: The only admitting subgroups are ['T_semidirect_{1}_cyclic_translation', 'T_semidirect_{1,11}_D12'].
- `killing_subgroups`: The subgroups ['T_semidirect_{1,5}_unit5_axis', 'T_semidirect_{1,7}_unit7_axis', 'T_semidirect_{1,5,7,11}_full_affine_Aut'] all mix opposite branch signs in one quotient orbit.
- `unit_obstruction`: The obstruction appears exactly when the admitted unit subgroup contains unit 5 or unit 7; unit 11 alone is harmless because it preserves the chi_11 sign pairs.
- `full_aut_consequence`: Full affine Aut contains both unit 5 and unit 7, so any full-Aut invariant support source is constant on the mixed branch orbit and cannot export chi_11 polarity.
- `strict_limit`: This subgroup-lattice fact does not derive the missing unit-axis premise and does not discharge QW-2191.

## Hard limits

- No identity K_legacy_ont == K_strict_gate is asserted.
- No legacy physical-role transfer onto K_strict_gate is used.
- No theorem derives chi_11, shell-labels, unit-axis bit, exact-cover clauses, or cardinality 5 from strict geometry.
- No QW-2191 discharge is claimed.
- No ToE closure is claimed.
- Result classifies affine subgroup quotients; it does not supply an internal strict source for selecting the D12/unit-axis premise.
