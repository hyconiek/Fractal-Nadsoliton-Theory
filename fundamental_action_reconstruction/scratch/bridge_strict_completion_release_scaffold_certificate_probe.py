#!/usr/bin/env python3
"""Scratch probe: release scaffold certificate for strict-kernel bridge documents.

This probe is deliberately editorial *and* computational.  It checks that the
new release-scaffold Markdown files exist, contain the required bridge,
Lagrangian/EOM, and role-transfer guardrails, and agree with the current finite
certificate ledger.  It does not claim a new bridge theorem.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FAR = ROOT / "fundamental_action_reconstruction"
OUT_JSON = HERE / "bridge_strict_completion_release_scaffold_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_release_scaffold_certificate_report.md"

SCAFFOLD_DOCS = {
    "strict_kernel_transformation_diagrams": FAR / "DIAGRAMS_STRICT_KERNEL_TRANSFORMATION.md",
    "strict_kernel_lagrangian_eom_draft": FAR / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
    "strict_kernel_role_transfer_audit_draft": FAR / "STRICT_KERNEL_LEGACY_ROLE_TRANSFER_AUDIT_DRAFT.md",
}

SOURCE_REPORTS = {
    "chain_integrity": HERE / "bridge_strict_completion_certificate_chain_integrity_report.json",
    "finite_bridge_assembly": HERE / "bridge_strict_completion_legacy_to_strict_finite_bridge_assembly_certificate_report.json",
    "symbolic_cancellation": HERE / "bridge_strict_completion_legacy_to_strict_symbolic_cancellation_certificate_report.json",
    "role_transfer_pre_audit": HERE / "bridge_strict_completion_legacy_role_transfer_pre_audit_certificate_report.json",
    "role_transfer_minimal_obligation_lattice": HERE / "bridge_strict_completion_legacy_role_transfer_minimal_obligation_lattice_certificate_report.json",
    "theorem_frontier_low_weight_extension": HERE / "bridge_strict_completion_theorem_frontier_low_weight_extension_certificate_report.json",
    "anchor_h1_generator_classification": HERE / "bridge_strict_completion_anchor_h1_generator_classification_certificate_report.json",
    "strict_lagrangian_symbolic_export_p1866": FAR / "generated/p1866_s816_strict_kernel_to_full_lagrangian_and_eom_symbolic_export.json",
    "strict_schematic_eom_spectrum_p2315": FAR / "generated/p2315_s1265_strict_schematic_lagrangian_eom_kernel_spectrum_probe.json",
}

REQUIRED_SNIPPETS = {
    "strict_kernel_transformation_diagrams": [
        "K_legacy_ont(d)",
        "K_strict_gate(d)",
        "K_legacy_ont(d) * A(d) * P(d) * D(d) = K_strict_gate(d)",
        "nonlinear damping/compression",
        "GF(2) phase pattern",
        "nadsoliton -> light -> matter -> emergent observer",
        "No unqualified identity",
        "No legacy physical-role transfer",
        "No `QW-2191` selector discharge",
        "No ToE closure",
    ],
    "strict_kernel_lagrangian_eom_draft": [
        "K_strict_gate(d) = cos(omega*d + phi)/(1 + beta*d^eta)",
        "beta = 1",
        "eta  = 9/5",
        "c0 = K_strict(1)",
        "m2_eff = m2*(1 + c0)",
        "L_total = L_scalar + L_fermion + L_gauge + L_gravity",
        "Box Phi + m2_eff*Phi + lam_eff*Phi^3",
        "eom_phi_proxy_1d",
        "QW-2191",
        "does not close QW-2191 or ToE",
    ],
    "strict_kernel_role_transfer_audit_draft": [
        "sin^2(theta_W)=alpha_geo/12",
        "alpha_EM^-1 = alpha_geo/(2*beta_tors)*(1-beta_tors)",
        "beta^N",
        "beta_tors -> chi_11",
        "nonlinear damping/compression",
        "GF(2) phase-sign solution",
        "all four audited roles",
        "No legacy role is transferred",
        "No QW-2191 selector discharge",
        "No ToE closure",
    ],
}


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing source report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def read_doc(path: Path) -> str:
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8")


def contains_all(text: str, snippets: list[str]) -> bool:
    return all(snippet in text for snippet in snippets)


def build_payload() -> dict[str, Any]:
    reports = {name: load_json(path) for name, path in SOURCE_REPORTS.items()}
    docs = {name: read_doc(path) for name, path in SCAFFOLD_DOCS.items()}

    doc_checks = {
        name: {
            "path": rel(SCAFFOLD_DOCS[name]),
            "exists": SCAFFOLD_DOCS[name].exists(),
            "required_snippet_count": len(REQUIRED_SNIPPETS[name]),
            "matched_required_snippet_count": sum(snippet in docs[name] for snippet in REQUIRED_SNIPPETS[name]),
            "all_required_snippets_present": contains_all(docs[name], REQUIRED_SNIPPETS[name]),
            "line_count": len(docs[name].splitlines()),
        }
        for name in SCAFFOLD_DOCS
    }

    chain_summary = reports["chain_integrity"]["chain_summary"]
    finite_bridge_summary = reports["finite_bridge_assembly"]["finite_bridge_assembly_summary"]
    symbolic_summary = reports["symbolic_cancellation"]["symbolic_cancellation_summary"]
    role_summary = reports["role_transfer_pre_audit"]["role_transfer_pre_audit_summary"]
    lattice_summary = reports["role_transfer_minimal_obligation_lattice"]["role_transfer_minimal_obligation_summary"]
    low_weight_summary = reports["theorem_frontier_low_weight_extension"]["theorem_frontier_low_weight_extension_summary"]
    anchor_summary = reports["anchor_h1_generator_classification"]["classification_summary"]
    p1866 = reports["strict_lagrangian_symbolic_export_p1866"]
    p2315 = reports["strict_schematic_eom_spectrum_p2315"]

    source_consistency_checks = {
        "chain_integrity_component_witnesses_present": chain_summary["legacy_to_strict_finite_bridge_assembly_certified"] and chain_summary["legacy_role_transfer_pre_audit_certified"] and chain_summary["anchor_h1_generator_classification_certified"],
        "finite_bridge_reconstructs_without_role_transfer": finite_bridge_summary["assembled_map_reconstructs_strict_kernel"] and finite_bridge_summary["comparison_scope_only"] and not finite_bridge_summary["legacy_role_transfer_exported"],
        "symbolic_cancellation_formula_only": symbolic_summary["symbolic_cancellation_formula_exported"] and not symbolic_summary["strict_dynamic_source_exported"] and not symbolic_summary["legacy_role_transfer_exported"],
        "role_transfer_pre_audit_blocks_all_roles": role_summary["all_roles_blocked_now"] and role_summary["roles_transferred_now"] == 0 and not role_summary["role_transfer_theorem_exported"],
        "role_transfer_lattice_has_four_atom_obligation": lattice_summary["global_minimal_obligation_size"] == 4 and not lattice_summary["role_transfer_theorem_exported"],
        "low_weight_frontier_keeps_bridge_role_toe_open": low_weight_summary["no_singleton_closes_bridge_role_or_toe"] and low_weight_summary["no_pair_closes_bridge"] and low_weight_summary["no_pair_closes_role_transfer"] and low_weight_summary["no_pair_closes_toe"],
        "anchor_h1_type_audit_preserves_selector_gap": anchor_summary["left_anchor_is_c0_gauge_fix_not_c1_generator"] and anchor_summary["selector_source_remains_open"],
        "p1866_symbolic_lagrangian_export_loaded": p1866["status"] == "OPEN_OBSTRUCTION_WITH_TRACE" and "L_total_decomposition" in p1866["full_lagrangian_non_skeleton"],
        "p2315_schematic_eom_keeps_qw2191_open": p2315["gatekeeper_checks"]["schematic_eom_derived"] and p2315["gatekeeper_checks"]["no_qw2191_discharge_claimed"],
        "chain_already_certifies_component_witnesses": chain_summary["legacy_to_strict_finite_bridge_assembly_certified"] and chain_summary["legacy_role_transfer_pre_audit_certified"] and chain_summary["anchor_h1_generator_classification_certified"],
    }

    cross_checks = {
        "three_scaffold_files_present": all(check["exists"] for check in doc_checks.values()),
        "strict_diagram_has_bridge_history": doc_checks["strict_kernel_transformation_diagrams"]["all_required_snippets_present"],
        "lagrangian_eom_draft_has_p1866_exports": doc_checks["strict_kernel_lagrangian_eom_draft"]["all_required_snippets_present"],
        "role_transfer_audit_blocks_all_legacy_roles": doc_checks["strict_kernel_role_transfer_audit_draft"]["all_required_snippets_present"],
        **source_consistency_checks,
        "no_identity_bridge_claimed": "No unqualified identity" in docs["strict_kernel_transformation_diagrams"] and "K_legacy_ont == K_strict_gate" in docs["strict_kernel_transformation_diagrams"],
        "no_role_transfer_theorem_claimed": "No legacy role is transferred" in docs["strict_kernel_role_transfer_audit_draft"] and role_summary["roles_transferred_now"] == 0,
        "no_qw2191_discharge": "No QW-2191 selector discharge" in docs["strict_kernel_role_transfer_audit_draft"] and p2315["gatekeeper_checks"]["no_qw2191_discharge_claimed"],
        "no_toe_closure": "No ToE closure" in docs["strict_kernel_role_transfer_audit_draft"] and not role_summary["toe_closure_claimed"],
    }

    all_cross_checks_pass = all(cross_checks.values())

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_RELEASE_SCAFFOLD_CERTIFICATE__THREE_MD_FILES_AUDITED_NO_CLOSURE",
        "status": "release-scaffold-files-present-and-consistent-with-finite-bridge-ledger-no-false-pass",
        "source_reports": {name: rel(path) for name, path in SOURCE_REPORTS.items()},
        "scaffold_docs": {name: rel(path) for name, path in SCAFFOLD_DOCS.items()},
        "doc_checks": doc_checks,
        "source_consistency_checks": source_consistency_checks,
        "cross_checks": cross_checks,
        "all_cross_checks_pass": all_cross_checks_pass,
        "release_scaffold_summary": {
            "strict_kernel_diagram_scaffold_ready": cross_checks["strict_diagram_has_bridge_history"],
            "strict_lagrangian_eom_scaffold_ready": cross_checks["lagrangian_eom_draft_has_p1866_exports"],
            "strict_role_transfer_audit_scaffold_ready": cross_checks["role_transfer_audit_blocks_all_legacy_roles"],
            "finite_bridge_ledger_inherited": source_consistency_checks["finite_bridge_reconstructs_without_role_transfer"],
            "symbolic_lagrangian_eom_exports_inherited": source_consistency_checks["p1866_symbolic_lagrangian_export_loaded"] and source_consistency_checks["p2315_schematic_eom_keeps_qw2191_open"],
            "all_legacy_roles_blocked_pending_theorems": source_consistency_checks["role_transfer_pre_audit_blocks_all_roles"],
            "no_identity_bridge_claimed": cross_checks["no_identity_bridge_claimed"],
            "no_role_transfer_theorem_claimed": cross_checks["no_role_transfer_theorem_claimed"],
            "no_qw2191_discharge": cross_checks["no_qw2191_discharge"],
            "no_toe_closure": cross_checks["no_toe_closure"],
        },
        "proof_certificate": {
            "strict_diagram_step": "DIAGRAMS_STRICT_KERNEL_TRANSFORMATION.md records K_legacy_ont as intermediate, K_strict_gate as enriched strict kernel, and the explicit finite APD bridge equation while preserving no raw identity, no role-transfer, no QW-2191, and no ToE limits.",
            "lagrangian_eom_step": "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md imports the P1866 K_strict -> (c0,c1,c2) -> effective couplings -> L_total scaffold and P2315 schematic EOM/spectrum status, but keeps 4D residual-zero EOM and selector closure open.",
            "role_transfer_step": "STRICT_KERNEL_LEGACY_ROLE_TRANSFER_AUDIT_DRAFT.md classifies sin^2(theta_W), alpha_EM^-1, beta^N, and beta_tors->chi_11 as blocked pending alpha_geo, beta_tors, beta-power, and chi11 theorem atoms; no legacy role is transferred now.",
            "ledger_step": "The scaffold is backed by the existing finite bridge assembly, symbolic cancellation, role-transfer pre-audit/lattice, anchor/H1 audit, and theorem-frontier reports, so this is a release-build package rather than a new theorem closure.",
        },
        "hard_limits": [
            "No unqualified identity K_legacy_ont == K_strict_gate is claimed.",
            "No strict dynamical derivation of A/P/D, omega/phi, beta/eta, or chi_11 is exported.",
            "No legacy physical-role transfer theorem is claimed.",
            "No beta_tors -> chi_11 theorem is claimed.",
            "No QW-2191 selector discharge is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# Strict-kernel release scaffold certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Scaffold docs",
        "",
    ]
    for name, check in payload["doc_checks"].items():
        lines.append(f"- `{name}`: exists=`{check['exists']}`, snippets=`{check['matched_required_snippet_count']}/{check['required_snippet_count']}`, lines=`{check['line_count']}`")
    lines.extend(["", "## Cross-checks", ""])
    for key, value in payload["cross_checks"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Summary", ""])
    for key, value in payload["release_scaffold_summary"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Hard limits", ""])
    for limit in payload["hard_limits"]:
        lines.append(f"- {limit}")
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    payload = build_payload()
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_markdown(payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
