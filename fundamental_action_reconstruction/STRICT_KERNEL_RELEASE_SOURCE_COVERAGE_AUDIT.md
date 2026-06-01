# Strict kernel release source coverage audit

Status: `DRAFT_V0_REPO_GREP_PROVENANCE_AND_SOURCE_COVERAGE__NO_THEOREM_CLOSURE`

This file records the repo-grep provenance for the strict-kernel release scaffold. It exists to prevent accidental duplication and to make clear which older strict/legacy artifacts are being folded into the new release-facing Markdown files.

## 1. Grep strategy used before extending the scaffold

The working search pass was intentionally broad and bridge-focused:

```text
rg -n "DIAGRAMS_KERNEL_TRANSFORMATION|DIAGRAMS_STRICT_KERNEL_TRANSFORMATION|STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT|STRICT_KERNEL_LEGACY_ROLE_TRANSFER_AUDIT_DRAFT|P1622|P1866|P2315|p2316|strict.*Lagrangian|full strict Lagrangian|EOM|legacy.*strict|strict.*legacy|role-transfer|role transfer|sin\^2\(theta_W\)|alpha_EM|beta_tors.*chi|chi_11|QW-2191" fundamental_action_reconstruction DIAGRAMS_KERNEL_TRANSFORMATION.md RELEASE_8_NB_CLOSURE_AUDIT_EN_PL.md --glob '!*.json' | head -500
```

This search found four source families that matter for the next release: legacy kernel history, strict Lagrangian/EOM exports, legacy-to-strict bridge certificates, and role-transfer/selector boundary theorems.

## 2. Source-family coverage map

| Family | Representative source files | How the new scaffold uses them |
|---|---|---|
| Legacy kernel history | `DIAGRAMS_KERNEL_TRANSFORMATION.md`; `K1_LEGACY_ONTOLOGICAL_KERNEL_VS_STRICT_GATE_KERNEL_SPLIT_NOTE.md`; `S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md` | Establishes `K_legacy_ont` as restored intermediate bridge kernel, not a discarded or final strict kernel. |
| Strict kernel and full Lagrangian/EOM route | `P1622_S572_FULL_STRICT_LAGRANGIAN_DENSITY_AND_EOM_PACKET_PL.md`; `p1866_s816_strict_kernel_to_full_lagrangian_and_eom_symbolic_export.py`; `generated/p1866_s816_strict_kernel_to_full_lagrangian_and_eom_symbolic_export.json`; `p2315_s1265_strict_schematic_lagrangian_eom_kernel_spectrum_probe.py`; `p2316_s1266_strict_current_best_lagrangian_eom_coverage_audit_probe.py`; `generated/p2316_s1266_strict_current_best_lagrangian_eom_coverage_audit_probe.json` | Establishes that a working strict Lagrangian/EOM route exists, but full tensor-resolved theorem-grade EOM closure remains open. |
| Legacy -> strict finite bridge ledger | `bridge_strict_completion_legacy_to_strict_finite_bridge_assembly_certificate_probe.py`; `bridge_strict_completion_legacy_to_strict_symbolic_cancellation_certificate_probe.py`; `bridge_strict_completion_certificate_chain_integrity_probe.py`; `bridge_strict_completion_release_scaffold_certificate_probe.py` | Establishes the APD/diagonal finite comparison ledger, not a strict dynamical bridge theorem. |
| Role-transfer and selector boundaries | `STRICT_KERNEL_LEGACY_ROLE_TRANSFER_AUDIT_DRAFT.md`; `N103_CURRENT_STRICT_SIDE_GRAVITY_HIERARCHY_ROLE_EQUIVALENCE_BOUNDARY_THEOREM.md`; `n87_current_strict_side_fine_structure_role_equivalence_boundary_theorem.py`; `bridge_strict_completion_legacy_role_transfer_pre_audit_certificate_probe.py`; `bridge_strict_completion_legacy_role_transfer_minimal_obligation_lattice_certificate_probe.py`; `bridge_strict_completion_anchor_h1_generator_classification_certificate_probe.py`; `bridge_strict_completion_theorem_frontier_low_weight_extension_certificate_probe.py` | Establishes that strict-side candidate objects may exist while legacy semantic role-transfer verdicts remain missing. |

## 3. Non-duplication decision

The new release-facing files are not meant to replace the source probes. Their job is narrower:

1. `DIAGRAMS_STRICT_KERNEL_TRANSFORMATION.md` is the human release diagram scaffold for the bridge story.
2. `STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md` is the human release scaffold for the strict Lagrangian/EOM route.
3. `STRICT_KERNEL_LEGACY_ROLE_TRANSFER_AUDIT_DRAFT.md` is the human release scaffold for updating legacy physical roles.
4. This file records grep/source provenance so later tex/md release writing can cite the underlying probe families rather than rediscovering them.

## 4. Coverage findings

The grep/source audit supports the following honest coverage statement:

- The legacy history side is covered enough for a draft bridge narrative.
- The finite APD bridge ledger is covered enough for a draft strict transformation diagram.
- The strict Lagrangian/EOM side has a working route and symbolic exports, but the P2316 audit still says the full tensor-resolved theorem-grade EOM/Lagrangian needed for Task-3 is missing.
- The role-transfer side has candidate strict objects and pre-audits, but no legacy semantic transfer theorem is exported.
- The selector side still keeps `QW-2191` open; the anchor/H1 audit only disambiguates where the bit lives and where it does not live.

## 5. Known release gaps after this grep pass

Before turning the scaffolds into final release tex/md, the next drafts must fold in more explicit references to:

1. the P1622 sector-level strict Lagrangian packet;
2. the P1866 symbolic coefficient and EOM-proxy export;
3. the P2315 schematic EOM/spectrum non-selector result;
4. the P2316 current-best strict Lagrangian/EOM coverage audit;
5. the strict-side fine-structure and gravity-hierarchy boundary theorems (`N87`/`N103` families);
6. the finite bridge assembly and symbolic cancellation reports;
7. the role-transfer pre-audit and minimal obligation lattice;
8. the theorem-frontier low-weight result showing that `chi11_selector_source` is the only singleton unlock but still not a bridge/role/ToE closure.

## 6. Hard limits

- This source-coverage audit does not prove a bridge theorem.
- This source-coverage audit does not prove full tensor-resolved EOM closure.
- This source-coverage audit does not transfer legacy physical roles.
- This source-coverage audit does not prove `beta_tors -> chi_11`.
- This source-coverage audit does not discharge `QW-2191`.
- This source-coverage audit does not close ToE.

## 7. Exact required source-path ledger for the certificate probe

The automated source-coverage certificate checks the following exact paths:

```text
DIAGRAMS_KERNEL_TRANSFORMATION.md
fundamental_action_reconstruction/K1_LEGACY_ONTOLOGICAL_KERNEL_VS_STRICT_GATE_KERNEL_SPLIT_NOTE.md
fundamental_action_reconstruction/S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md
fundamental_action_reconstruction/P1622_S572_FULL_STRICT_LAGRANGIAN_DENSITY_AND_EOM_PACKET_PL.md
fundamental_action_reconstruction/p1866_s816_strict_kernel_to_full_lagrangian_and_eom_symbolic_export.py
fundamental_action_reconstruction/generated/p1866_s816_strict_kernel_to_full_lagrangian_and_eom_symbolic_export.json
fundamental_action_reconstruction/p2315_s1265_strict_schematic_lagrangian_eom_kernel_spectrum_probe.py
fundamental_action_reconstruction/p2316_s1266_strict_current_best_lagrangian_eom_coverage_audit_probe.py
fundamental_action_reconstruction/generated/p2316_s1266_strict_current_best_lagrangian_eom_coverage_audit_probe.json
fundamental_action_reconstruction/scratch/bridge_strict_completion_legacy_to_strict_finite_bridge_assembly_certificate_probe.py
fundamental_action_reconstruction/scratch/bridge_strict_completion_legacy_to_strict_symbolic_cancellation_certificate_probe.py
fundamental_action_reconstruction/scratch/bridge_strict_completion_certificate_chain_integrity_probe.py
fundamental_action_reconstruction/scratch/bridge_strict_completion_release_scaffold_certificate_probe.py
fundamental_action_reconstruction/STRICT_KERNEL_LEGACY_ROLE_TRANSFER_AUDIT_DRAFT.md
fundamental_action_reconstruction/N103_CURRENT_STRICT_SIDE_GRAVITY_HIERARCHY_ROLE_EQUIVALENCE_BOUNDARY_THEOREM.md
fundamental_action_reconstruction/n87_current_strict_side_fine_structure_role_equivalence_boundary_theorem.py
fundamental_action_reconstruction/scratch/bridge_strict_completion_legacy_role_transfer_pre_audit_certificate_probe.py
fundamental_action_reconstruction/scratch/bridge_strict_completion_legacy_role_transfer_minimal_obligation_lattice_certificate_probe.py
fundamental_action_reconstruction/scratch/bridge_strict_completion_anchor_h1_generator_classification_certificate_probe.py
fundamental_action_reconstruction/scratch/bridge_strict_completion_theorem_frontier_low_weight_extension_certificate_probe.py
```

This exact-path ledger is part of the repo-grep provenance certificate; it is a source map, not a theorem export.
