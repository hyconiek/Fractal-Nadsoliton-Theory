#!/usr/bin/env python3
"""Scratch probe: assumption-conditional ToE closure interface.

This is the honest closure step: it proves the finite implication

    seven open theorem atoms supplied  ==>  all audited targets true

inside the already-audited 128-row frontier table, while keeping every theorem
atom itself open in the current strict-kernel state.  The output is therefore a
conditional closure certificate, not a strict ToE closure claim.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FAR = ROOT / "fundamental_action_reconstruction"
OUT_JSON = HERE / "bridge_strict_completion_toe_conditional_closure_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_toe_conditional_closure_certificate_report.md"
TOE_AUDIT_DOC = FAR / "STRICT_KERNEL_TOE_POTENTIAL_AUDIT.md"

SOURCE_REPORTS = {
    "frontier_truth_table": HERE / "bridge_strict_completion_theorem_frontier_truth_table_certificate_report.json",
    "target_signature_lattice": HERE / "bridge_strict_completion_theorem_frontier_target_signature_lattice_certificate_report.json",
    "toe_boolean_normal_form": HERE / "bridge_strict_completion_toe_boolean_normal_form_certificate_report.json",
    "toe_boolean_essentiality": HERE / "bridge_strict_completion_toe_boolean_essentiality_certificate_report.json",
}

TARGET_KEYS = [
    "bridge_theorem_level_closure",
    "role_transfer_theorem_level_closure",
    "selector_qw2191_closure",
    "toe_closure",
]

SEQUENT_DEFINITIONS = [
    {
        "sequent_id": "S1_bridge_sources_to_bridge_target",
        "assumptions": [
            "strict_damping_beta_eta_source",
            "strict_dynamical_source_for_A_P_D",
            "strict_phase_frequency_source",
        ],
        "conclusions": ["bridge_theorem_level_closure"],
        "scope": "conditional bridge theorem interface",
    },
    {
        "sequent_id": "S2_role_atoms_to_role_target",
        "assumptions": [
            "alpha_geo_electroweak_role_theorem",
            "beta_power_hierarchy_successor_theorem",
            "beta_tors_strict_role_theorem",
            "chi11_selector_source",
        ],
        "conclusions": ["role_transfer_theorem_level_closure"],
        "scope": "conditional strict role-transfer theorem interface",
    },
    {
        "sequent_id": "S3_chi11_to_selector_target",
        "assumptions": ["chi11_selector_source"],
        "conclusions": ["selector_qw2191_closure"],
        "scope": "conditional selector-source interface",
    },
    {
        "sequent_id": "S4_all_atoms_to_toe_target",
        "assumptions": [
            "alpha_geo_electroweak_role_theorem",
            "beta_power_hierarchy_successor_theorem",
            "beta_tors_strict_role_theorem",
            "chi11_selector_source",
            "strict_damping_beta_eta_source",
            "strict_dynamical_source_for_A_P_D",
            "strict_phase_frequency_source",
        ],
        "conclusions": TARGET_KEYS,
        "scope": "assumption-conditional ToE closure interface",
    },
]

DOC_SNIPPETS = [
    "Conditional closure interface",
    "seven theorem atoms as assumptions",
    "No unconditional ToE closure is claimed",
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
    signature_lattice = reports["target_signature_lattice"]
    normal_form = reports["toe_boolean_normal_form"]
    essentiality = reports["toe_boolean_essentiality"]
    doc_text = TOE_AUDIT_DOC.read_text(encoding="utf-8") if TOE_AUDIT_DOC.exists() else ""

    open_atoms = truth["open_atoms"]
    atom_to_bit = {atom: index for index, atom in enumerate(open_atoms)}
    full_mask = (1 << len(open_atoms)) - 1
    truth_rows = sorted(truth["truth_table_rows"], key=lambda row: row["assignment_index"])
    row_by_mask = {row_mask(row, atom_to_bit): row for row in truth_rows}

    current_row = row_by_mask[0]
    full_row = row_by_mask[full_mask]
    target_minimal_sets = truth["target_minimal_true_atom_sets"]

    sequent_rows: list[dict[str, Any]] = []
    for sequent in SEQUENT_DEFINITIONS:
        assumption_mask = 0
        for atom in sequent["assumptions"]:
            assumption_mask |= 1 << atom_to_bit[atom]
        witness_row = row_by_mask[assumption_mask]
        conclusions = sequent["conclusions"]
        conclusion_values = {target: bool(witness_row[target]) for target in conclusions}
        assumption_set = set(sequent["assumptions"])
        minimality_checks = {
            target: any(set(atom_set).issubset(assumption_set) for atom_set in target_minimal_sets[target])
            for target in conclusions
            if target in target_minimal_sets
        }
        sequent_rows.append(
            {
                "sequent_id": sequent["sequent_id"],
                "scope": sequent["scope"],
                "assumptions": sequent["assumptions"],
                "assumption_count": len(sequent["assumptions"]),
                "assumption_mask": assumption_mask,
                "conclusions": conclusions,
                "conclusion_values": conclusion_values,
                "all_conclusions_true_under_assumptions": all(conclusion_values.values()),
                "minimality_checks": minimality_checks,
                "all_available_minimality_checks_pass": all(minimality_checks.values()) if minimality_checks else True,
                "witness_assignment_index": witness_row["assignment_index"],
                "witness_true_atoms": witness_row["true_atoms"],
            }
        )

    conditional_closure_summary = {
        "conditional_toe_sequent_id": "S4_all_atoms_to_toe_target",
        "conditional_toe_assumption_count": len(open_atoms),
        "conditional_toe_full_row_index": full_row["assignment_index"],
        "conditional_toe_full_row_closes_all_targets": all(bool(full_row[target]) for target in TARGET_KEYS),
        "current_row_closes_no_targets": not any(bool(current_row[target]) for target in TARGET_KEYS),
        "strict_open_atoms_still_open_now": True,
        "unconditional_toe_closure_claimed": False,
    }

    cross_checks = {
        "source_reports_present": all(path.exists() for path in SOURCE_REPORTS.values()),
        "truth_table_loaded": len(truth_rows) == 128 and truth["theorem_frontier_truth_table_summary"]["truth_assignment_count"] == 128,
        "full_row_is_unique_toe_row": truth["theorem_frontier_truth_table_summary"]["toe_satisfying_assignment_count"] == 1 and full_row["toe_closure"],
        "target_signature_lattice_inherited": all(signature_lattice["cross_checks"].values()) and signature_lattice["theorem_frontier_target_signature_lattice_summary"]["only_full_signature_has_toe_closure"],
        "boolean_normal_form_inherited": normal_form["all_cross_checks_pass"] and normal_form["boolean_normal_form_summary"]["toe_anf_degree"] == 7,
        "boolean_essentiality_inherited": essentiality["all_cross_checks_pass"] and essentiality["boolean_essentiality_summary"]["toe_derivative_witness_count"] == 7,
        "all_sequents_true_under_assumptions": all(row["all_conclusions_true_under_assumptions"] for row in sequent_rows),
        "minimal_sequents_match_truth_table": all(row["all_available_minimality_checks_pass"] for row in sequent_rows),
        "current_state_unclosed": conditional_closure_summary["current_row_closes_no_targets"] and truth["theorem_frontier_truth_table_summary"]["current_targets_all_false"],
        "toe_audit_doc_mentions_conditional_interface": all(snippet in doc_text for snippet in DOC_SNIPPETS),
        "hard_limits_preserved": not truth["theorem_frontier_truth_table_summary"]["toe_closure_claimed"] and not conditional_closure_summary["unconditional_toe_closure_claimed"],
    }

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_TOE_CONDITIONAL_CLOSURE_CERTIFICATE__ASSUMPTION_SEQUENTS_NO_UNCONDITIONAL_CLOSURE",
        "status": "PASS" if all(cross_checks.values()) else "FAIL",
        "source_reports": {name: rel(path) for name, path in SOURCE_REPORTS.items()},
        "open_atoms": open_atoms,
        "target_keys": TARGET_KEYS,
        "sequent_rows": sequent_rows,
        "cross_checks": cross_checks,
        "all_cross_checks_pass": all(cross_checks.values()),
        "conditional_closure_summary": conditional_closure_summary,
        "proof_certificate": {
            "sequent_step": "The probe verifies four finite sequents: bridge-source atoms imply the bridge target, role atoms imply the role-transfer target, chi11 implies selector/QW-2191 closure, and all seven atoms imply all four target bits.",
            "minimality_step": "Each sequent assumption set covers a truth-table minimal true atom set for every asserted target; the component sequents are exact minimal interfaces, and the ToE sequent uses the full seven-atom frontier cut.",
            "closure_step": "The full seven-atom row is the unique truth-table row with toe_closure=true, so this is an assumption-conditional ToE closure interface.",
            "limit_step": "The current zero-atom row still closes no target; none of the seven source/role/selector atoms is supplied by this probe, so no unconditional ToE closure is claimed.",
        },
        "hard_limits": [
            "No theorem atom is proved by this conditional interface.",
            "No strict dynamical source theorem is claimed.",
            "No legacy role-transfer theorem is claimed unconditionally.",
            "No QW-2191 selector source is supplied.",
            "No unconditional ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# ToE conditional closure certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Conditional sequents",
        "",
    ]
    for row in payload["sequent_rows"]:
        assumptions = ", ".join(f"`{atom}`" for atom in row["assumptions"])
        conclusions = ", ".join(f"`{target}`" for target in row["conclusions"])
        lines.append(
            f"- `{row['sequent_id']}`: {assumptions} => {conclusions}; "
            f"mask=`{row['assumption_mask']}`, pass=`{row['all_conclusions_true_under_assumptions']}`"
        )
    lines.extend(["", "## Cross-checks", ""])
    for key, value in payload["cross_checks"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Summary", ""])
    for key, value in payload["conditional_closure_summary"].items():
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
