# P3140/S2090 Axiom-augmented selector premise calculus

Status: `P3140_AXIOM_AUGMENTED_SELECTOR_PREMISE_CALCULUS_NON_STRICT`

## Constructed object
- `S_sel^ax(D_HL) axiom-augmented selector premise calculus`
- Axioms:
  - `A_origin`: postulate one absolute support representative r0 in Z12
  - `A_lambda`: postulate one polarity lambda0 in {+1,-1}
  - `A_coupling`: postulate that the selected pair (r0,lambda0) couples to the P3134 D_HL formula as the selector channel

## Finite theorem
`P3140_T1_axiom_augmented_DHL_selector_calculus`: On the finite pair space Z12 x {±1}, no current strict artifact selects a unique joint pair.  Adding A_origin and A_lambda reduces the pair space to one stipulated pair, and adding A_coupling turns that pair into a conditional D_HL selector channel.  This is a valid non-strict axiom-augmented closure, but it is not a strict-source theorem and does not discharge QW-2191 in the strict core.

## Finite counts
- `pair_space_size`: `24`
- `no_axiom_orbit_size_under_translation_aut_lambda_flip`: `24`
- `axiom_packages_tested`: `8`
- `minimal_non_strict_closing_packages`: `1`

## Obligations
- `absolute_origin`: strict_current=`False`, axiom_needed=`A_origin`, non_strict_if_assumed=`True`
- `unpaired_lambda`: strict_current=`False`, axiom_needed=`A_lambda`, non_strict_if_assumed=`True`
- `DHL_selector_coupling`: strict_current=`False`, axiom_needed=`A_coupling`, non_strict_if_assumed=`True`
- `strict_source_provenance`: strict_current=`False`, axiom_needed=`None`, non_strict_if_assumed=`False`
- `variational_unit_Ltotal_bridge`: strict_current=`False`, axiom_needed=`downstream theorem even after selector axioms`, non_strict_if_assumed=`False`

## Decision
Assuming axioms can honestly break the selector symmetry, but only as a labelled non-strict extension.  The minimal D_HL package is A_origin + A_lambda + A_coupling: it selects one stipulated pair and couples it to D_HL.  It does not prove strict provenance, does not erase QW-2191 from the strict core, and does not by itself supply variational/unit/L_total or ToE closure.

## Recommendation
If the project accepts an axiom route, formalize a separate non-strict axiom-augmented branch with A_origin, A_lambda, and A_coupling as explicit assumptions, then audit downstream consequences without calling them strict.  If the project stays strict-core, do not use the axiom package as evidence; pivot to one genuinely new strict source object that can derive one of these axioms rather than assume it.
