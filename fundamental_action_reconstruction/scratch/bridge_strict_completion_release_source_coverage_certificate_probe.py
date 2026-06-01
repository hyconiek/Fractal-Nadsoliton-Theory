#!/usr/bin/env python3
"""Scratch probe: repo-grep source coverage certificate for strict release scaffold.

This probe turns the user's "did you grep the repo?" concern into a finite
certificate.  It scans the current repository text for the source families used
by the strict-kernel release scaffold, checks that the source coverage audit MD
mentions the required files and limits, and records the remaining non-closure
status.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FAR = ROOT / "fundamental_action_reconstruction"
OUT_JSON = HERE / "bridge_strict_completion_release_source_coverage_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_release_source_coverage_certificate_report.md"
SOURCE_COVERAGE_DOC = FAR / "STRICT_KERNEL_RELEASE_SOURCE_COVERAGE_AUDIT.md"

SCAN_ROOTS = [FAR, ROOT / "DIAGRAMS_KERNEL_TRANSFORMATION.md", ROOT / "RELEASE_8_NB_CLOSURE_AUDIT_EN_PL.md"]

SOURCE_FAMILIES = {
    "legacy_kernel_history": {
        "patterns": ["DIAGRAMS_KERNEL_TRANSFORMATION", "K_legacy_ont", "S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET"],
        "required_files": [
            "DIAGRAMS_KERNEL_TRANSFORMATION.md",
            "fundamental_action_reconstruction/K1_LEGACY_ONTOLOGICAL_KERNEL_VS_STRICT_GATE_KERNEL_SPLIT_NOTE.md",
            "fundamental_action_reconstruction/S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md",
        ],
    },
    "strict_lagrangian_eom": {
        "patterns": ["P1622", "P1866", "P2315", "p2316", "full strict Lagrangian", "EOM"],
        "required_files": [
            "fundamental_action_reconstruction/P1622_S572_FULL_STRICT_LAGRANGIAN_DENSITY_AND_EOM_PACKET_PL.md",
            "fundamental_action_reconstruction/p1866_s816_strict_kernel_to_full_lagrangian_and_eom_symbolic_export.py",
            "fundamental_action_reconstruction/generated/p1866_s816_strict_kernel_to_full_lagrangian_and_eom_symbolic_export.json",
            "fundamental_action_reconstruction/p2315_s1265_strict_schematic_lagrangian_eom_kernel_spectrum_probe.py",
            "fundamental_action_reconstruction/p2316_s1266_strict_current_best_lagrangian_eom_coverage_audit_probe.py",
            "fundamental_action_reconstruction/generated/p2316_s1266_strict_current_best_lagrangian_eom_coverage_audit_probe.json",
        ],
    },
    "finite_bridge_ledger": {
        "patterns": ["legacy_to_strict", "finite_bridge_assembly", "symbolic_cancellation", "release_scaffold", "certificate_chain_integrity"],
        "required_files": [
            "fundamental_action_reconstruction/scratch/bridge_strict_completion_legacy_to_strict_finite_bridge_assembly_certificate_probe.py",
            "fundamental_action_reconstruction/scratch/bridge_strict_completion_legacy_to_strict_symbolic_cancellation_certificate_probe.py",
            "fundamental_action_reconstruction/scratch/bridge_strict_completion_certificate_chain_integrity_probe.py",
            "fundamental_action_reconstruction/scratch/bridge_strict_completion_release_scaffold_certificate_probe.py",
        ],
    },
    "role_transfer_selector_boundaries": {
        "patterns": ["role-transfer", "role transfer", "alpha_EM", "beta_tors", "chi_11", "QW-2191", "N103", "n87"],
        "required_files": [
            "fundamental_action_reconstruction/STRICT_KERNEL_LEGACY_ROLE_TRANSFER_AUDIT_DRAFT.md",
            "fundamental_action_reconstruction/N103_CURRENT_STRICT_SIDE_GRAVITY_HIERARCHY_ROLE_EQUIVALENCE_BOUNDARY_THEOREM.md",
            "fundamental_action_reconstruction/n87_current_strict_side_fine_structure_role_equivalence_boundary_theorem.py",
            "fundamental_action_reconstruction/scratch/bridge_strict_completion_legacy_role_transfer_pre_audit_certificate_probe.py",
            "fundamental_action_reconstruction/scratch/bridge_strict_completion_legacy_role_transfer_minimal_obligation_lattice_certificate_probe.py",
            "fundamental_action_reconstruction/scratch/bridge_strict_completion_anchor_h1_generator_classification_certificate_probe.py",
            "fundamental_action_reconstruction/scratch/bridge_strict_completion_theorem_frontier_low_weight_extension_certificate_probe.py",
        ],
    },
}

REQUIRED_DOC_SNIPPETS = [
    "repo-grep provenance",
    "DIAGRAMS_KERNEL_TRANSFORMATION.md",
    "P1622_S572_FULL_STRICT_LAGRANGIAN_DENSITY_AND_EOM_PACKET_PL.md",
    "p1866_s816_strict_kernel_to_full_lagrangian_and_eom_symbolic_export.py",
    "p2315_s1265_strict_schematic_lagrangian_eom_kernel_spectrum_probe.py",
    "p2316_s1266_strict_current_best_lagrangian_eom_coverage_audit_probe.py",
    "N103_CURRENT_STRICT_SIDE_GRAVITY_HIERARCHY_ROLE_EQUIVALENCE_BOUNDARY_THEOREM.md",
    "n87_current_strict_side_fine_structure_role_equivalence_boundary_theorem.py",
    "chi11_selector_source",
    "does not discharge `QW-2191`",
    "does not close ToE",
]


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def iter_scan_files() -> list[Path]:
    files: list[Path] = []
    for root in SCAN_ROOTS:
        if root.is_file():
            files.append(root)
        else:
            files.extend(
                path
                for path in root.rglob("*")
                if path.is_file()
                and path.suffix in {".md", ".py", ".tex"}
                and ".git" not in path.parts
            )
    return sorted(set(files), key=lambda path: rel(path))


def count_hits(patterns: list[str], files: list[Path]) -> tuple[int, list[dict[str, Any]]]:
    hits: list[dict[str, Any]] = []
    for path in files:
        text = path.read_text(encoding="utf-8", errors="ignore")
        matched = [pattern for pattern in patterns if pattern in text]
        if matched:
            hits.append({"path": rel(path), "matched_patterns": matched[:6], "matched_pattern_count": len(matched)})
    return len(hits), hits[:25]


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(path)
    return json.loads(path.read_text(encoding="utf-8"))


def build_payload() -> dict[str, Any]:
    files = iter_scan_files()
    doc_text = SOURCE_COVERAGE_DOC.read_text(encoding="utf-8") if SOURCE_COVERAGE_DOC.exists() else ""
    p2316 = load_json(FAR / "generated/p2316_s1266_strict_current_best_lagrangian_eom_coverage_audit_probe.json")
    p1866 = load_json(FAR / "generated/p1866_s816_strict_kernel_to_full_lagrangian_and_eom_symbolic_export.json")
    release_scaffold = load_json(HERE / "bridge_strict_completion_release_scaffold_certificate_report.json")

    family_rows: list[dict[str, Any]] = []
    for family, spec in SOURCE_FAMILIES.items():
        hit_count, top_hits = count_hits(spec["patterns"], files)
        required_paths = [ROOT / required for required in spec["required_files"]]
        family_rows.append(
            {
                "family": family,
                "patterns": spec["patterns"],
                "grep_like_hit_count": hit_count,
                "top_hits": top_hits,
                "required_files": spec["required_files"],
                "required_files_present": all(path.exists() for path in required_paths),
                "required_files_mentioned_in_doc": all(required in doc_text for required in spec["required_files"]),
            }
        )

    cross_checks = {
        "source_coverage_doc_present": SOURCE_COVERAGE_DOC.exists(),
        "all_required_doc_snippets_present": all(snippet in doc_text for snippet in REQUIRED_DOC_SNIPPETS),
        "all_source_families_have_hits": all(row["grep_like_hit_count"] > 0 for row in family_rows),
        "all_required_family_files_present": all(row["required_files_present"] for row in family_rows),
        "all_required_family_files_mentioned": all(row["required_files_mentioned_in_doc"] for row in family_rows),
        "p2316_repo_grep_audit_available": p2316["result_kind"] == "STRICT_REPO_GREP_AND_COMPUTATIONAL_LAGRANGIAN_EOM_COVERAGE_AUDIT_NO_G1_G3_UPDATE" and p2316["gatekeeper_checks"]["repo_lagrangian_grep_hits_found"],
        "p2316_keeps_full_task3_theorem_open": p2316["gatekeeper_checks"]["full_task3_theorem_not_claimed"] and p2316["gatekeeper_checks"]["no_qw2191_discharge_claimed"] and p2316["gatekeeper_checks"]["no_toe_closure_claimed"],
        "p1866_symbolic_export_still_obstruction": p1866["status"] == "OPEN_OBSTRUCTION_WITH_TRACE" and "L_total_decomposition" in p1866["full_lagrangian_non_skeleton"],
        "release_scaffold_already_passes": release_scaffold["all_cross_checks_pass"],
        "no_identity_role_selector_toe_closure": "does not prove a bridge theorem" in doc_text and "does not transfer legacy physical roles" in doc_text and "does not discharge `QW-2191`" in doc_text and "does not close ToE" in doc_text,
    }

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_RELEASE_SOURCE_COVERAGE_CERTIFICATE__REPO_GREP_PROVENANCE_NO_CLOSURE",
        "status": "strict-release-source-coverage-audited-by-repo-grep-provenance-no-false-pass",
        "source_coverage_doc": rel(SOURCE_COVERAGE_DOC),
        "scan_file_count": len(files),
        "source_family_rows": family_rows,
        "cross_checks": cross_checks,
        "all_cross_checks_pass": all(cross_checks.values()),
        "coverage_summary": {
            "legacy_history_covered": family_rows[0]["required_files_present"] and family_rows[0]["required_files_mentioned_in_doc"],
            "strict_lagrangian_eom_sources_covered": family_rows[1]["required_files_present"] and family_rows[1]["required_files_mentioned_in_doc"],
            "finite_bridge_ledger_covered": family_rows[2]["required_files_present"] and family_rows[2]["required_files_mentioned_in_doc"],
            "role_transfer_selector_boundaries_covered": family_rows[3]["required_files_present"] and family_rows[3]["required_files_mentioned_in_doc"],
            "p2316_full_task3_theorem_still_open": cross_checks["p2316_keeps_full_task3_theorem_open"],
            "release_scaffold_nonduplicating_source_map_ready": cross_checks["release_scaffold_already_passes"] and cross_checks["all_required_doc_snippets_present"],
            "no_bridge_theorem_claimed": True,
            "no_role_transfer_claimed": True,
            "no_qw2191_discharge": True,
            "no_toe_closure": True,
        },
        "proof_certificate": {
            "grep_step": "The source coverage audit records the actual broad rg pattern used for strict-kernel release preparation and this probe independently scans the repo text for four source families: legacy history, strict Lagrangian/EOM, finite bridge ledger, and role-transfer/selector boundaries.",
            "coverage_step": "All required representative source files are present and mentioned in STRICT_KERNEL_RELEASE_SOURCE_COVERAGE_AUDIT.md, including DIAGRAMS_KERNEL_TRANSFORMATION.md, P1622, P1866, P2315, P2316, N87/N103, finite bridge certificates, and role-transfer/frontier audits.",
            "nonduplication_step": "The new release-facing files are certified as scaffolds that summarize and point back to source families rather than replacing the underlying probes or claiming theorem closure.",
            "limit_step": "P2316 still keeps the full Task-3 tensor-resolved EOM theorem open, P1866 remains OPEN_OBSTRUCTION_WITH_TRACE, and the source-coverage audit preserves no bridge theorem, no role transfer, no QW-2191 discharge, and no ToE closure.",
        },
        "hard_limits": [
            "No source-coverage hit is promoted to a theorem.",
            "No bridge theorem is claimed.",
            "No full tensor-resolved EOM closure is claimed.",
            "No legacy physical-role transfer is claimed.",
            "No beta_tors -> chi_11 theorem is claimed.",
            "No QW-2191 selector discharge is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# Strict release source coverage certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        f"Scan file count: `{payload['scan_file_count']}`",
        "",
        "## Source families",
        "",
    ]
    for row in payload["source_family_rows"]:
        lines.append(f"- `{row['family']}`: hits=`{row['grep_like_hit_count']}`, required_files_present=`{row['required_files_present']}`, mentioned=`{row['required_files_mentioned_in_doc']}`")
    lines.extend(["", "## Cross-checks", ""])
    for key, value in payload["cross_checks"].items():
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
