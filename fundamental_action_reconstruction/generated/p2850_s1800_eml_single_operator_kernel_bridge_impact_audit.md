# P2850/S1800 EML single-operator kernel bridge impact audit

Status: `P2850_EML_SINGLE_OPERATOR_KERNEL_BRIDGE_IMPACT_AUDIT_NO_CLOSURE`

## External sources checked
- HTML: https://arxiv.org/html/2603.21852v2
- PDF: https://arxiv.org/pdf/2603.21852
- relevant_claim_summary=EML(x,y)=Exp[x]-Log[y], with terminal 1, is claimed to reconstruct the elementary scientific-calculator basis; the paper frames this as syntactic/operator-basis universality and separates heuristic search from independent verification.

## EML sanity check
- max_exp_identity_abs_error=0.0

## Elementary syntax classification
- legacy_formula_elementary=True
- strict_formula_elementary=True
- eml_representation_relevance=syntax_only_unless_typed_source_law_added

## Kernel value separation
- max_legacy_strict_kernel_abs_gap_on_audited_distances=2.6915403157237625

## Premise matrix
- changes_p2849_bridge_boundary=False
- missing_premises=['eml_exports_parameter_source_law', 'eml_exports_beta_eta_source', 'eml_exports_amplitude_source', 'eml_exports_phase_topological_selector', 'eml_exports_completion_map_semantics', 'eml_exports_role_transfer_theorem']

## Boundary
The arXiv EML result changes the available expression syntax, not the typed source obligations.  It supports treating legacy and strict kernel formulas as elementary expressions that may be encoded in one operator basis, but it does not export beta/eta, amplitude, phase/topological selector, completion-map semantics, unit-bearing L_total coupling, or role-transfer source laws.  P2849's damping/compression obstruction is unchanged.

## Recommendation
Use the EML paper only as an optional syntax/compression lens.  The next proof-grade move remains one typed bridge-source atom: either a new strict source law for eta and target-independent beta, or a separate amplitude-normalization passage theorem.  If no such typed source premise is supplied, preserve the no-new-live-frontier certificate.
