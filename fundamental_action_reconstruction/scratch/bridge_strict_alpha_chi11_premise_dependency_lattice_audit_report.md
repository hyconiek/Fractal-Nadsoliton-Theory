# Strict-alpha chi_11 premise-dependency lattice audit

- Result kind: `SCRATCH_STRICT_ALPHA_CHI11_PREMISE_DEPENDENCY_LATTICE_AUDIT__NO_STRICT_FULL_AUT_SOURCE`
- Status: `minimal-premise-antichain-computed-for-existing-chi11-support-audits`
- Premise subsets checked: `256`.
- Input reports loaded: `7`.
- Full-Aut exports chi_11 polarity: `False`.
- D12 chi_11 rank: `13`.
- Cyclic chi_11 dimension: `13`.

## Minimal premise antichain

| Outcome | Minimal premise count | Minimal premise sets |
|---|---:|---|
| `locate_branch_full_aut_block_without_polarity` | `1` | `[['full_aut_unoriented_block_amplitude']]` |
| `host_nonzero_cyclic_chi11_character_space` | `2` | `[['cyclic_translation_quotient', 'chi11_unit_character_choice']]` |
| `host_d12_chi11_covariant_module` | `2` | `[['d12_reduced_quotient', 'chi11_unit_character_choice']]` |
| `branch_normalized_d12_chi11_family` | `3` | `[['d12_reduced_quotient', 'chi11_unit_character_choice', 'branch_normalization_constraint']]` |
| `unique_branch_generator_by_sparsest_extension` | `4` | `[['d12_reduced_quotient', 'chi11_unit_character_choice', 'branch_normalization_constraint', 'sparsest_extension_selector']]` |
| `unique_branch_generator_by_max_shell_imbalance` | `5` | `[['d12_reduced_quotient', 'chi11_unit_character_choice', 'branch_normalization_constraint', 'shell_label_d1_d5_axis', 'max_shell_imbalance_selector']]` |
| `strict_full_aut_internal_chi11_polarity_source` | `None` | `[]` |

## Proof certificate

- This audit reuses the already-generated C(12,5)=792 support reports and checks their shared finite-model counts before building a premise antichain.
- The affine subgroup lattice says full Aut does not admit chi_11; the full-Aut block-amplitude report locates the branch block but has exports_chi11_polarity=false; the nonhistogram audit has zero full-Aut singleton A5-not-A1 classifiers.
- A nonzero chi_11 carrier appears only after quotient/premise reduction: cyclic translation quotient plus chi_11 character gives dimension 13, and D12 plus chi_11 character gives the 13 two-cycle module.
- The branch-normalized D12 chi_11 family is not unique until an extra selector is imported.  Current finite audits certify two conditional selectors: sparsest extension and shell-labelled max imbalance.
- No row in this premise lattice supplies an internal strict full-Aut source for the chi_11 polarity, the shell-label d1/d5 axis, or QW-2191 selector closure.

## Hard limits

- No identity K_legacy_ont == K_strict_gate is asserted.
- No legacy physical role is transferred onto K_strict_gate.
- No theorem derives chi_11 from full-Aut strict geometry.
- No theorem derives the shell-label d1/d5 axis or unit-axis bit from strict full-Aut support data.
- No theorem derives the sparsest-extension selector or max-shell-imbalance selector as strict-core closure.
- No QW-2191 discharge is claimed.
- No ToE closure is claimed.
