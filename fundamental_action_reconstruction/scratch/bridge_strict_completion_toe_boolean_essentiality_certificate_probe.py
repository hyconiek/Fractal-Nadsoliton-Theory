#!/usr/bin/env python3
"""Scratch probe: Boolean essentiality/derivative audit for ToE closure.

This probe computes finite Boolean derivatives for the audited theorem-frontier
targets.  For ToE, every open atom is essential exactly at the six-atom nearest
miss where all other atoms are true.  Thus each atom has a concrete witness for
necessity, without exporting any theorem-level closure.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FAR = ROOT / "fundamental_action_reconstruction"
OUT_JSON = HERE / "bridge_strict_completion_toe_boolean_essentiality_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_toe_boolean_essentiality_certificate_report.md"
TOE_AUDIT_DOC = FAR / "STRICT_KERNEL_TOE_POTENTIAL_AUDIT.md"

SOURCE_REPORTS = {
    "frontier_truth_table": HERE / "bridge_strict_completion_theorem_frontier_truth_table_certificate_report.json",
    "toe_boolean_normal_form": HERE / "bridge_strict_completion_toe_boolean_normal_form_certificate_report.json",
    "toe_proper_subset_obstruction": HERE / "bridge_strict_completion_toe_proper_subset_obstruction_certificate_report.json",
}

TARGET_KEYS = [
    "bridge_theorem_level_closure",
    "role_transfer_theorem_level_closure",
    "selector_qw2191_closure",
    "toe_closure",
]

EXPECTED_DERIVATIVE_SUPPORT_COUNTS = {
    "bridge_theorem_level_closure": {
        "strict_dynamical_source_for_A_P_D": 16,
        "strict_phase_frequency_source": 16,
        "strict_damping_beta_eta_source": 16,
    },
    "role_transfer_theorem_level_closure": {
        "alpha_geo_electroweak_role_theorem": 8,
        "beta_power_hierarchy_successor_theorem": 8,
        "beta_tors_strict_role_theorem": 8,
        "chi11_selector_source": 8,
    },
    "selector_qw2191_closure": {
        "chi11_selector_source": 64,
    },
    "toe_closure": {
        "alpha_geo_electroweak_role_theorem": 1,
        "beta_power_hierarchy_successor_theorem": 1,
        "beta_tors_strict_role_theorem": 1,
        "chi11_selector_source": 1,
        "strict_damping_beta_eta_source": 1,
        "strict_dynamical_source_for_A_P_D": 1,
        "strict_phase_frequency_source": 1,
    },
}

DOC_SNIPPETS = [
    "Boolean essentiality",
    "Boolean derivative",
    "No ToE closure is claimed",
]


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(path)
    return json.loads(path.read_text(encoding="utf-8"))


def row_mask(row: dict[str, Any], atom_to_bit: dict[str, int]) -> int:
    mask = 0
    for atom in row["true_atoms"]:
        mask |= 1 << atom_to_bit[atom]
    return mask


def build_payload() -> dict[str, Any]:
    reports = {name: load_json(path) for name, path in SOURCE_REPORTS.items()}
    truth = reports["frontier_truth_table"]
    normal_form = reports["toe_boolean_normal_form"]
    proper_subset = reports["toe_proper_subset_obstruction"]
    doc_text = TOE_AUDIT_DOC.read_text(encoding="utf-8") if TOE_AUDIT_DOC.exists() else ""

    open_atoms = truth["open_atoms"]
    atom_to_bit = {atom: index for index, atom in enumerate(open_atoms)}
    truth_rows = sorted(truth["truth_table_rows"], key=lambda row: row["assignment_index"])
    row_by_mask = {row_mask(row, atom_to_bit): row for row in truth_rows}
    full_mask = (1 << len(open_atoms)) - 1

    derivative_rows: list[dict[str, Any]] = []
    for target in TARGET_KEYS:
        expected_for_target = EXPECTED_DERIVATIVE_SUPPORT_COUNTS[target]
        for atom in open_atoms:
            bit = 1 << atom_to_bit[atom]
            support_masks: list[int] = []
            examples: list[dict[str, Any]] = []
            for base_mask in range(1 << len(open_atoms)):
                if base_mask & bit:
                    continue
                row0 = row_by_mask[base_mask]
                row1 = row_by_mask[base_mask | bit]
                if bool(row0[target]) ^ bool(row1[target]):
                    support_masks.append(base_mask)
                    if len(examples) < 3:
                        examples.append(
                            {
                                "base_mask": base_mask,
                                "raised_mask": base_mask | bit,
                                "base_true_atoms": row0["true_atoms"],
                                "raised_true_atoms": row1["true_atoms"],
                                "base_value": row0[target],
                                "raised_value": row1[target],
                            }
                        )
            derivative_rows.append(
                {
                    "target": target,
                    "atom": atom,
                    "support_count": len(support_masks),
                    "expected_support_count": expected_for_target.get(atom, 0),
                    "is_essential_for_target": len(support_masks) > 0,
                    "support_masks": support_masks,
                    "examples": examples,
                }
            )

    toe_derivative_rows = [row for row in derivative_rows if row["target"] == "toe_closure"]
    nearest_miss_by_atom = {
        row["missing_atom"]: row for row in proper_subset["nearest_miss_rows"]
    }
    toe_nearest_miss_rows = []
    for row in toe_derivative_rows:
        atom = row["atom"]
        base_mask = row["support_masks"][0] if row["support_masks"] else None
        toe_nearest_miss_rows.append(
            {
                "atom": atom,
                "support_count": row["support_count"],
                "base_mask": base_mask,
                "expected_base_mask": full_mask ^ (1 << atom_to_bit[atom]),
                "nearest_miss_signature": nearest_miss_by_atom[atom]["target_signature"],
                "nearest_miss_explanation": nearest_miss_by_atom[atom]["deficit_explanation"],
            }
        )

    cross_checks = {
        "source_reports_present": all(path.exists() for path in SOURCE_REPORTS.values()),
        "truth_table_loaded": len(truth_rows) == 128 and truth["theorem_frontier_truth_table_summary"]["truth_assignment_count"] == 128,
        "normal_form_inherited": normal_form["all_cross_checks_pass"] and normal_form["boolean_normal_form_summary"]["toe_anf_degree"] == 7,
        "proper_subset_obstruction_inherited": proper_subset["all_cross_checks_pass"] and proper_subset["proper_subset_obstruction_summary"]["all_127_proper_subsets_fail_toe"],
        "all_derivative_counts_match_expected": all(row["support_count"] == row["expected_support_count"] for row in derivative_rows),
        "toe_each_atom_essential_once": all(row["support_count"] == 1 for row in toe_derivative_rows),
        "toe_derivative_supports_are_nearest_misses": all(row["base_mask"] == row["expected_base_mask"] for row in toe_nearest_miss_rows),
        "nonparticipants_have_zero_derivative": all(
            row["support_count"] == 0
            for row in derivative_rows
            if row["atom"] not in EXPECTED_DERIVATIVE_SUPPORT_COUNTS[row["target"]]
        ),
        "toe_audit_doc_mentions_essentiality": all(snippet in doc_text for snippet in DOC_SNIPPETS),
        "hard_limits_preserved": not truth["theorem_frontier_truth_table_summary"]["toe_closure_claimed"] and not normal_form["boolean_normal_form_summary"]["toe_closure_claimed"] and not proper_subset["proper_subset_obstruction_summary"]["toe_closure_claimed"],
    }

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_TOE_BOOLEAN_ESSENTIALITY_CERTIFICATE__SEVEN_DERIVATIVE_WITNESSES_NO_CLOSURE",
        "status": "PASS" if all(cross_checks.values()) else "FAIL",
        "source_reports": {name: rel(path) for name, path in SOURCE_REPORTS.items()},
        "open_atoms": open_atoms,
        "target_keys": TARGET_KEYS,
        "derivative_rows": derivative_rows,
        "toe_nearest_miss_rows": toe_nearest_miss_rows,
        "cross_checks": cross_checks,
        "all_cross_checks_pass": all(cross_checks.values()),
        "boolean_essentiality_summary": {
            "toe_derivative_witness_count": sum(row["support_count"] for row in toe_derivative_rows),
            "toe_each_atom_essential_once": cross_checks["toe_each_atom_essential_once"],
            "toe_derivative_supports_are_nearest_misses": cross_checks["toe_derivative_supports_are_nearest_misses"],
            "component_derivative_counts_match_expected": cross_checks["all_derivative_counts_match_expected"],
            "nonparticipant_derivatives_zero": cross_checks["nonparticipants_have_zero_derivative"],
            "toe_closure_claimed": False,
        },
        "proof_certificate": {
            "derivative_step": "The finite Boolean derivative audit computes f(x with atom=0) XOR f(x with atom=1) for every target and open atom across the 128-row truth table.",
            "toe_step": "For toe_closure, each of the seven atoms has derivative support count 1, and that support is exactly the six-atom nearest miss where all other atoms are true.",
            "component_step": "The component derivative support counts match the monomial degrees: bridge participants have 16 supports, role-transfer participants 8, selector chi11 64, and all nonparticipants 0.",
            "limit_step": "This is an essentiality witness over the audited Boolean frontier; it exports no source theorem, no bridge theorem, no role-transfer theorem, no QW-2191 discharge, and no ToE closure.",
        },
        "hard_limits": [
            "No Boolean derivative witness is promoted to a source theorem.",
            "No bridge theorem is claimed.",
            "No role-transfer theorem is claimed.",
            "No QW-2191 selector discharge is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# ToE Boolean essentiality certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## ToE nearest-miss derivative witnesses",
        "",
    ]
    for row in payload["toe_nearest_miss_rows"]:
        lines.append(
            f"- `{row['atom']}`: support_count=`{row['support_count']}`, "
            f"base_mask=`{row['base_mask']}`, signature=`{row['nearest_miss_signature']}`"
        )
    lines.extend(["", "## Cross-checks", ""])
    for key, value in payload["cross_checks"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Summary", ""])
    for key, value in payload["boolean_essentiality_summary"].items():
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
