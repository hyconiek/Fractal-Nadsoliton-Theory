# Strict alpha chi_11 candidate-source transform audit

Status: `candidate-records-do-not-export-chi11-without-shell-label-or-reduced-symmetry-import`

## Candidate summary

- Candidate count: `7`
- Full-Aut invariant candidate count: `3`
- Pair-distinguishing candidate count: `4`
- Numeric chi_11-covariant candidate count: `2`
- Allowed strict chi_11 source candidate count: `0`
- chi_11-covariant but imported candidates: `['shell_labelled_d5_minus_d1_count', 'shell_labelled_energy_difference_(d5+d6)-(d1+d6)']`

## Proof certificate

- `finite_domain`: The audit is over the four unit branches A1,A5,A7,A11 and the exact unit action of Aut(Z_12)={1,5,7,11}.
- `required_transform`: A chi_11 source must be constant on reversal pairs, separate contiguous_pair from d5_pair, and flip under units 5 and 7 with kernel {1,11}.
- `invariant_no_go`: Full-Aut invariant candidates that distinguish the branch pairs: []; hence invariant scalar records do not provide the bit.
- `covariant_imports`: Candidates with chi_11 covariance are ['shell_labelled_d5_minus_d1_count', 'shell_labelled_energy_difference_(d5+d6)-(d1+d6)'], but each requires a d1-vs-d5 shell-label or reduced-symmetry premise.
- `allowed_source_count`: Allowed strict chi_11 source candidates found without imported shell labels: 0.
- `boundary`: The computation classifies candidate records; it does not prove that the listed candidate family exhausts all possible strict nadsoliton sources.

## Candidate rows

- `support_size`: kind=`full_aut_invariant_scalar`, distinguishes_pairs=`False`, chi11_covariant=`False`, allowed_source=`False`
- `full_aut_orbit_histogram_O15_O2_O3_O4_O6`: kind=`full_aut_invariant_scalar`, distinguishes_pairs=`False`, chi11_covariant=`False`, allowed_source=`False`
- `antipodal_pair_count_d6`: kind=`full_aut_invariant_scalar`, distinguishes_pairs=`False`, chi11_covariant=`False`, allowed_source=`False`
- `dihedral_gap_necklace`: kind=`dihedral_invariant_not_full_aut`, distinguishes_pairs=`True`, chi11_covariant=`False`, allowed_source=`False`
- `raw_distance_histogram_d1_to_d6`: kind=`shell_labelled_record`, distinguishes_pairs=`True`, chi11_covariant=`False`, allowed_source=`False`
- `shell_labelled_d5_minus_d1_count`: kind=`chi11_covariant_but_imports_shell_label`, distinguishes_pairs=`True`, chi11_covariant=`True`, allowed_source=`False`
- `shell_labelled_energy_difference_(d5+d6)-(d1+d6)`: kind=`chi11_covariant_but_imports_shell_label`, distinguishes_pairs=`True`, chi11_covariant=`True`, allowed_source=`False`

## Hard limits

- No identity K_legacy_ont == K_strict_gate is used or claimed.
- No legacy role transfer to K_strict_gate, alpha_geo, beta_tors, or D_f is made.
- No theorem derives the chi_11-kernel, shell-label d1 vs d5, unit-axis bit, exact-cover clauses, or cardinality 5 from strict nadsoliton geometry.
- The result is a finite candidate-source transform audit over listed branch records, not an exhaustive strict-source theorem.
- No QW-2191 discharge.
- No ToE closure.
