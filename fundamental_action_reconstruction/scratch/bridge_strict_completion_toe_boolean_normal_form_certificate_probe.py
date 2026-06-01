#!/usr/bin/env python3
"""Scratch probe: Boolean normal-form audit for strict-kernel ToE closure.

This probe computes the GF(2) algebraic normal form (ANF/Zhegalkin form) of
all frontier target bits and checks that the ToE target is the single degree-7
monomial containing every open theorem atom.  It is a finite normal-form audit,
not a theorem export.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FAR = ROOT / "fundamental_action_reconstruction"
OUT_JSON = HERE / "bridge_strict_completion_toe_boolean_normal_form_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_toe_boolean_normal_form_certificate_report.md"
TOE_AUDIT_DOC = FAR / "STRICT_KERNEL_TOE_POTENTIAL_AUDIT.md"

SOURCE_REPORTS = {
    "frontier_truth_table": HERE / "bridge_strict_completion_theorem_frontier_truth_table_certificate_report.json",
    "toe_proper_subset_obstruction": HERE / "bridge_strict_completion_toe_proper_subset_obstruction_certificate_report.json",
    "toe_potential_readiness": HERE / "bridge_strict_completion_toe_potential_readiness_certificate_report.json",
}

TARGET_KEYS = [
    "bridge_theorem_level_closure",
    "role_transfer_theorem_level_closure",
    "selector_qw2191_closure",
    "toe_closure",
]

EXPECTED_DEGREES = {
    "bridge_theorem_level_closure": 3,
    "role_transfer_theorem_level_closure": 4,
    "selector_qw2191_closure": 1,
    "toe_closure": 7,
}

DOC_SNIPPETS = [
    "Boolean normal form",
    "degree-7 monomial",
    "No ToE closure is claimed",
]


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(path)
    return json.loads(path.read_text(encoding="utf-8"))


def mobius_anf(values: list[int], variable_count: int) -> list[int]:
    coeffs = values[:]
    for bit in range(variable_count):
        step = 1 << bit
        for mask in range(1 << variable_count):
            if mask & step:
                coeffs[mask] ^= coeffs[mask ^ step]
    return coeffs


def monomial_atoms(mask: int, row_by_index: dict[int, dict[str, Any]]) -> list[str]:
    return row_by_index[mask]["true_atoms"]


def build_payload() -> dict[str, Any]:
    reports = {name: load_json(path) for name, path in SOURCE_REPORTS.items()}
    truth = reports["frontier_truth_table"]
    proper_subset = reports["toe_proper_subset_obstruction"]
    toe_readiness = reports["toe_potential_readiness"]
    doc_text = TOE_AUDIT_DOC.read_text(encoding="utf-8") if TOE_AUDIT_DOC.exists() else ""

    open_atoms = truth["open_atoms"]
    variable_count = len(open_atoms)
    truth_rows = sorted(truth["truth_table_rows"], key=lambda row: row["assignment_index"])
    row_by_index = {row["assignment_index"]: row for row in truth_rows}
    full_mask = (1 << variable_count) - 1

    normal_form_rows: list[dict[str, Any]] = []
    for target in TARGET_KEYS:
        values = [1 if row[target] else 0 for row in truth_rows]
        coeffs = mobius_anf(values, variable_count)
        monomial_masks = [mask for mask, coeff in enumerate(coeffs) if coeff]
        monomials = [
            {
                "mask": mask,
                "degree": len(monomial_atoms(mask, row_by_index)),
                "atoms": monomial_atoms(mask, row_by_index),
            }
            for mask in monomial_masks
        ]
        normal_form_rows.append(
            {
                "target": target,
                "monomial_count": len(monomials),
                "max_degree": max((row["degree"] for row in monomials), default=0),
                "monomials": monomials,
                "has_lower_degree_terms": any(row["degree"] < EXPECTED_DEGREES[target] for row in monomials),
            }
        )

    rows_by_target = {row["target"]: row for row in normal_form_rows}
    toe_row = rows_by_target["toe_closure"]
    component_degrees = {row["target"]: row["max_degree"] for row in normal_form_rows}
    component_monomial_counts = {row["target"]: row["monomial_count"] for row in normal_form_rows}

    cross_checks = {
        "source_reports_present": all(path.exists() for path in SOURCE_REPORTS.values()),
        "truth_table_128_rows_loaded": len(truth_rows) == 128 and truth["theorem_frontier_truth_table_summary"]["truth_assignment_count"] == 128,
        "proper_subset_obstruction_inherited": proper_subset["all_cross_checks_pass"] and proper_subset["proper_subset_obstruction_summary"]["all_127_proper_subsets_fail_toe"],
        "toe_readiness_inherited": toe_readiness["all_cross_checks_pass"] and toe_readiness["toe_potential_readiness_summary"]["toe_requires_all_7_open_atoms"],
        "component_normal_forms_single_monomial": all(row["monomial_count"] == 1 for row in normal_form_rows),
        "component_degrees_match_expected": component_degrees == EXPECTED_DEGREES,
        "toe_anf_single_full_degree_monomial": toe_row["monomial_count"] == 1 and toe_row["monomials"][0]["mask"] == full_mask and toe_row["monomials"][0]["degree"] == variable_count,
        "toe_anf_has_no_lower_degree_terms": not toe_row["has_lower_degree_terms"],
        "toe_doc_mentions_boolean_normal_form": all(snippet in doc_text for snippet in DOC_SNIPPETS),
        "hard_limits_preserved": not truth["theorem_frontier_truth_table_summary"]["toe_closure_claimed"] and not proper_subset["proper_subset_obstruction_summary"]["toe_closure_claimed"] and not toe_readiness["toe_potential_readiness_summary"]["toe_closure_claimed"],
    }

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_TOE_BOOLEAN_NORMAL_FORM_CERTIFICATE__ANF_DEGREE_7_FULL_MONOMIAL_NO_CLOSURE",
        "status": "PASS" if all(cross_checks.values()) else "FAIL",
        "source_reports": {name: rel(path) for name, path in SOURCE_REPORTS.items()},
        "open_atoms": open_atoms,
        "variable_count": variable_count,
        "target_keys": TARGET_KEYS,
        "normal_form_rows": normal_form_rows,
        "component_degrees": component_degrees,
        "component_monomial_counts": component_monomial_counts,
        "cross_checks": cross_checks,
        "all_cross_checks_pass": all(cross_checks.values()),
        "boolean_normal_form_summary": {
            "toe_anf_degree": toe_row["max_degree"],
            "toe_anf_monomial_count": toe_row["monomial_count"],
            "toe_anf_atoms": toe_row["monomials"][0]["atoms"],
            "toe_has_lower_degree_terms": toe_row["has_lower_degree_terms"],
            "all_target_anfs_are_single_monomials": cross_checks["component_normal_forms_single_monomial"],
            "proper_subset_obstruction_inherited": cross_checks["proper_subset_obstruction_inherited"],
            "toe_closure_claimed": False,
        },
        "proof_certificate": {
            "anf_step": "The GF(2) Mobius transform of the ToE truth column has exactly one ANF term: the degree-7 monomial containing all seven open theorem atoms.",
            "component_step": "The component target ANFs are also single monomials: bridge has degree 3, role-transfer has degree 4, selector/QW-2191 has degree 1, and ToE has degree 7.",
            "lower_degree_step": "Because ToE has no lower-degree ANF terms, no lower-order algebraic combination of fewer frontier atoms is present in the audited Boolean target definition.",
            "limit_step": "This normal-form certificate is finite algebraic bookkeeping over the current truth table; it exports no source theorem, no bridge theorem, no role-transfer theorem, no selector discharge, and no ToE closure.",
        },
        "hard_limits": [
            "No Boolean normal-form term is promoted to a source theorem.",
            "No bridge theorem is claimed.",
            "No role-transfer theorem is claimed.",
            "No QW-2191 selector discharge is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# ToE Boolean normal-form certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Normal-form rows",
        "",
    ]
    for row in payload["normal_form_rows"]:
        lines.append(f"- `{row['target']}`: monomials=`{row['monomial_count']}`, max_degree=`{row['max_degree']}`")
        for monomial in row["monomials"]:
            lines.append(f"  - degree `{monomial['degree']}` atoms=`{monomial['atoms']}`")
    lines.extend(["", "## Cross-checks", ""])
    for key, value in payload["cross_checks"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Summary", ""])
    for key, value in payload["boolean_normal_form_summary"].items():
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
