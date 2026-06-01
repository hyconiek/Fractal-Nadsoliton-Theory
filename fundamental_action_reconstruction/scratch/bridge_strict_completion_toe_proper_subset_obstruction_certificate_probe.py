#!/usr/bin/env python3
"""Scratch probe: proper-subset obstruction for strict-kernel ToE closure.

This is a monotone Boolean audit over the already-exported theorem frontier.  It
checks every proper subset of the seven open atoms and proves, by finite
enumeration plus nearest-miss classification, that no six-or-fewer atom package
can honestly be promoted to ToE closure.
"""
from __future__ import annotations

import json
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FAR = ROOT / "fundamental_action_reconstruction"
OUT_JSON = HERE / "bridge_strict_completion_toe_proper_subset_obstruction_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_toe_proper_subset_obstruction_certificate_report.md"
TOE_AUDIT_DOC = FAR / "STRICT_KERNEL_TOE_POTENTIAL_AUDIT.md"

SOURCE_REPORTS = {
    "frontier_truth_table": HERE / "bridge_strict_completion_theorem_frontier_truth_table_certificate_report.json",
    "frontier_cut": HERE / "bridge_strict_completion_theorem_frontier_cut_certificate_report.json",
    "target_signature_lattice": HERE / "bridge_strict_completion_theorem_frontier_target_signature_lattice_certificate_report.json",
    "toe_potential_readiness": HERE / "bridge_strict_completion_toe_potential_readiness_certificate_report.json",
}

TARGET_KEYS = [
    "bridge_theorem_level_closure",
    "role_transfer_theorem_level_closure",
    "selector_qw2191_closure",
    "toe_closure",
]

DEFICIT_BY_ATOM = {
    "strict_dynamical_source_for_A_P_D": "bridge-source atom missing: APD source is absent, so bridge and ToE fail",
    "strict_phase_frequency_source": "bridge-source atom missing: phase/frequency source is absent, so bridge and ToE fail",
    "strict_damping_beta_eta_source": "bridge-source atom missing: damping beta/eta source is absent, so bridge and ToE fail",
    "chi11_selector_source": "selector and role atom missing: QW-2191 selector, role transfer, and ToE fail",
    "alpha_geo_electroweak_role_theorem": "role-transfer atom missing: electroweak role theorem is absent, so role transfer and ToE fail",
    "beta_tors_strict_role_theorem": "role-transfer atom missing: beta_tors strict successor theorem is absent, so role transfer and ToE fail",
    "beta_power_hierarchy_successor_theorem": "role-transfer atom missing: hierarchy successor theorem is absent, so role transfer and ToE fail",
}

DOC_SNIPPETS = [
    "all seven open atoms are jointly required",
    "No ToE closure is claimed",
]


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(path)
    return json.loads(path.read_text(encoding="utf-8"))


def signature(row: dict[str, Any]) -> str:
    return "".join("1" if row[key] else "0" for key in TARGET_KEYS)


def build_payload() -> dict[str, Any]:
    reports = {name: load_json(path) for name, path in SOURCE_REPORTS.items()}
    truth = reports["frontier_truth_table"]
    cut = reports["frontier_cut"]
    target_lattice = reports["target_signature_lattice"]
    toe_readiness = reports["toe_potential_readiness"]
    doc_text = TOE_AUDIT_DOC.read_text(encoding="utf-8") if TOE_AUDIT_DOC.exists() else ""

    open_atoms = truth["open_atoms"]
    full_atom_set = set(open_atoms)
    truth_rows = truth["truth_table_rows"]
    row_by_atoms = {frozenset(row["true_atoms"]): row for row in truth_rows}
    proper_rows = [row for row in truth_rows if row["true_atom_count"] < len(open_atoms)]
    full_row = row_by_atoms[frozenset(open_atoms)]

    proper_subset_count_by_weight = Counter(row["true_atom_count"] for row in proper_rows)
    max_proper_true_targets = max(sum(bool(row[key]) for key in TARGET_KEYS) for row in proper_rows)
    max_proper_rows = [row for row in proper_rows if sum(bool(row[key]) for key in TARGET_KEYS) == max_proper_true_targets]
    nearest_miss_rows = [row for row in proper_rows if row["true_atom_count"] == len(open_atoms) - 1]

    missing_atom_rows: list[dict[str, Any]] = []
    for atom in open_atoms:
        row = row_by_atoms[frozenset(full_atom_set - {atom})]
        missing_atom_rows.append(
            {
                "missing_atom": atom,
                "true_atom_count": row["true_atom_count"],
                "target_signature": signature(row),
                "bridge_theorem_level_closure": row["bridge_theorem_level_closure"],
                "role_transfer_theorem_level_closure": row["role_transfer_theorem_level_closure"],
                "selector_qw2191_closure": row["selector_qw2191_closure"],
                "toe_closure": row["toe_closure"],
                "deficit_explanation": DEFICIT_BY_ATOM[atom],
            }
        )

    signature_counts_by_weight: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for row in proper_rows:
        signature_counts_by_weight[str(row["true_atom_count"])][signature(row)] += 1
    signature_counts_by_weight = {
        weight: dict(sorted(counts.items())) for weight, counts in sorted(signature_counts_by_weight.items(), key=lambda item: int(item[0]))
    }

    proper_subset_count = len(proper_rows)
    expected_proper_subset_count = (2 ** len(open_atoms)) - 1
    cross_checks = {
        "source_reports_present": all(path.exists() for path in SOURCE_REPORTS.values()),
        "truth_table_count_pass": truth["theorem_frontier_truth_table_summary"]["truth_assignment_count"] == 128 and len(truth_rows) == 128,
        "proper_subset_count_pass": proper_subset_count == expected_proper_subset_count == 127,
        "full_set_only_toe_pass": full_row["toe_closure"] and all(not row["toe_closure"] for row in proper_rows),
        "nearest_miss_count_pass": len(nearest_miss_rows) == len(open_atoms) == 7,
        "nearest_miss_rows_fail_toe": all(not row["toe_closure"] for row in nearest_miss_rows),
        "missing_atom_rows_cover_frontier": sorted(row["missing_atom"] for row in missing_atom_rows) == sorted(open_atoms),
        "minimal_cut_matches_truth_table": sorted(cut["theorem_frontier_cut_summary"]["minimal_open_leaf_cut_for_toe"]) == sorted(open_atoms),
        "target_lattice_full_signature_only": target_lattice["theorem_frontier_target_signature_lattice_summary"]["only_full_signature_has_toe_closure"],
        "readiness_certificate_inherited": toe_readiness["all_cross_checks_pass"] and toe_readiness["toe_potential_readiness_summary"]["toe_requires_all_7_open_atoms"],
        "toe_audit_doc_mentions_joint_requirement": all(snippet in doc_text for snippet in DOC_SNIPPETS),
        "hard_limits_preserved": not truth["theorem_frontier_truth_table_summary"]["toe_closure_claimed"] and not toe_readiness["toe_potential_readiness_summary"]["toe_closure_claimed"],
    }

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_TOE_PROPER_SUBSET_OBSTRUCTION_CERTIFICATE__ALL_127_PROPER_SUBSETS_FAIL_NO_CLOSURE",
        "status": "PASS" if all(cross_checks.values()) else "FAIL",
        "source_reports": {name: rel(path) for name, path in SOURCE_REPORTS.items()},
        "open_atoms": open_atoms,
        "target_keys": TARGET_KEYS,
        "proper_subset_count": proper_subset_count,
        "proper_subset_count_by_weight": {str(k): proper_subset_count_by_weight[k] for k in sorted(proper_subset_count_by_weight)},
        "signature_counts_by_weight": signature_counts_by_weight,
        "max_proper_true_targets": max_proper_true_targets,
        "max_proper_target_signatures": sorted({signature(row) for row in max_proper_rows}),
        "nearest_miss_rows": missing_atom_rows,
        "cross_checks": cross_checks,
        "all_cross_checks_pass": all(cross_checks.values()),
        "proper_subset_obstruction_summary": {
            "all_127_proper_subsets_fail_toe": cross_checks["proper_subset_count_pass"] and cross_checks["full_set_only_toe_pass"],
            "nearest_miss_count": len(nearest_miss_rows),
            "six_atom_packages_fail_toe": cross_checks["nearest_miss_rows_fail_toe"],
            "all_seven_atoms_required": cross_checks["minimal_cut_matches_truth_table"] and cross_checks["readiness_certificate_inherited"],
            "max_proper_true_targets": max_proper_true_targets,
            "toe_closure_claimed": False,
        },
        "proof_certificate": {
            "enumeration_step": "The probe enumerates all 127 proper subsets of the seven open theorem atoms and checks each truth-table row; every proper subset has toe_closure=false.",
            "nearest_miss_step": "The seven six-atom nearest misses are classified by the single missing atom: missing bridge-source atoms disable bridge+ToE, missing role atoms disable role-transfer+ToE, and missing chi11 disables selector+role+ToE.",
            "signature_step": "No proper subset reaches ToE signature 1111; the largest proper signatures are partial bridge/role/selector states, so partial progress is not promoted to ToE closure.",
            "limit_step": "This is a finite obstruction/readiness certificate only: it proves a proper-subset no-go over the audited Boolean frontier and exports no new theorem atom, bridge theorem, role-transfer theorem, QW-2191 discharge, or ToE closure.",
        },
        "hard_limits": [
            "No proper subset is promoted to ToE closure.",
            "No new theorem atom is exported.",
            "No bridge theorem is claimed.",
            "No role-transfer theorem is claimed.",
            "No QW-2191 selector discharge is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# ToE proper-subset obstruction certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Summary",
        "",
    ]
    for key, value in payload["proper_subset_obstruction_summary"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Proper subset counts by weight", ""])
    for weight, count in payload["proper_subset_count_by_weight"].items():
        lines.append(f"- weight `{weight}`: `{count}` subsets")
    lines.extend(["", "## Six-atom nearest misses", ""])
    for row in payload["nearest_miss_rows"]:
        lines.append(f"- missing `{row['missing_atom']}` -> signature `{row['target_signature']}`; {row['deficit_explanation']}")
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
    print(json.dumps({"status": payload["status"], "all_cross_checks_pass": payload["all_cross_checks_pass"]}, sort_keys=True))


if __name__ == "__main__":
    main()
