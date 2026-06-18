# P2849/S1799 damping/compression kernel bridge-atom audit

Status: `P2849_DAMPING_COMPRESSION_KERNEL_BRIDGE_ATOM_AUDIT_NO_CLOSURE`

## Parameters
- legacy_beta_tors=1/100
- strict_beta=1
- strict_eta=9/5

## Two-point exact-match matrix
For two distinct positive distances a,b, exact equality beta*d^eta = beta_tors*d at both points forces eta=1 and beta=beta_tors; therefore the strict eta=9/5 nonlinear compression cannot be sourced by legacy linear damping alone.
- pair=[1, 2]: forced_eta=1; strict_eta_matches=False
- pair=[1, 3]: forced_eta=1; strict_eta_matches=False
- pair=[2, 3]: forced_eta=1; strict_eta_matches=False
- pair=[3, 6]: forced_eta=1; strict_eta_matches=False
- pair=[6, 12]: forced_eta=1; strict_eta_matches=False

## beta_eff stats for strict eta
- min=0.0013697931912643546
- max=0.01
- max_over_min=7.3003721027184705
- constant_across_distances=False

## Premise matrix
- accepted_as_damping_compression_bridge_atom=False
- missing_premises=['two_point_exact_completion_map', 'target_independent_positive_beta_source', 'eta_source_law', 'amplitude_normalization_compatibility', 'phase_frequency_compatibility', 'selector_or_topological_source_compatibility', 'role_transfer_theorem_available']

## Boundary
P2849 tests exactly one kernel bridge atom: legacy linear torsion damping to strict nonlinear beta/eta compression.  A two-point exact-match theorem forces eta=1 and beta=beta_tors for any exact legacy-linear source, while the strict value eta=9/5 fails every audited two-point row.  Matching legacy damping at strict eta would require beta_eff(d)=beta_tors*d^(1-eta), which varies across distances.  Current artifacts therefore do not export a target-independent beta source or eta source law.

## Recommendation
Do not claim kernel bridge, role transfer, L_total, EOM, Hamiltonian, or ToE closure from damping/compression similarity.  The next admissible move is either one new strict source law for eta and target-independent beta, or a different single kernel bridge atom such as amplitude normalization passage.  Without a new source premise, preserve the no-new-live-frontier certificate.
