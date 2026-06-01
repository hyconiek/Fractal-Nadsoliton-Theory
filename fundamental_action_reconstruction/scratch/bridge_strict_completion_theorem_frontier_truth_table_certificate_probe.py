#!/usr/bin/env python3
"""Scratch probe: theorem-frontier closure truth-table certificate.

The frontier-cut report identifies seven open theorem atoms.  This probe makes
the closure logic fully explicit by enumerating all 2^7 truth assignments for
those atoms and evaluating four target predicates: bridge theorem closure,
role-transfer closure, selector/QW-2191 closure, and ToE closure.

It is deliberately a readiness certificate, not a closure theorem: the current
assignment has all seven atoms false, so every theorem-level target remains
false.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_theorem_frontier_truth_table_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_theorem_frontier_truth_table_certificate_report.md"

SOURCE_REPORTS = {
    "theorem_frontier_cut": HERE / "bridge_strict_completion_theorem_frontier_cut_certificate_report.json",
    "role_transfer_minimal_obligation_lattice": HERE / "bridge_strict_completion_legacy_role_transfer_minimal_obligation_lattice_certificate_report.json",
    "component_gap_matrix": HERE / "bridge_strict_completion_legacy_to_strict_completion_component_gap_certificate_report.json",
    "anchor_h1_classification": HERE / "bridge_strict_completion_anchor_h1_generator_classification_certificate_report.json",
}

BRIDGE_ATOMS = [
    "strict_dynamical_source_for_A_P_D",
    "strict_phase_frequency_source",
    "strict_damping_beta_eta_source",
]
ROLE_ATOMS = [
    "alpha_geo_electroweak_role_theorem",
    "beta_tors_strict_role_theorem",
    "beta_power_hierarchy_successor_theorem",
    "chi11_selector_source",
]
SELECTOR_ATOMS = ["chi11_selector_source"]


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def all_assignments(atoms: list[str]) -> list[dict[str, bool]]:
    rows = []
    for bits in itertools.product([False, True], repeat=len(atoms)):
        rows.append(dict(zip(atoms, bits)))
    return rows


def target_values(assignment: dict[str, bool]) -> dict[str, bool]:
    bridge = all(assignment[atom] for atom in BRIDGE_ATOMS)
    role = all(assignment[atom] for atom in ROLE_ATOMS)
    selector = all(assignment[atom] for atom in SELECTOR_ATOMS)
    toe = bridge and role and selector
    return {
        "bridge_theorem_level_closure": bridge,
        "role_transfer_theorem_level_closure": role,
        "selector_qw2191_closure": selector,
        "toe_closure": toe,
    }


def true_atom_set(assignment: dict[str, bool], atoms: list[str]) -> list[str]:
    return [atom for atom in atoms if assignment[atom]]


def minimal_true_sets(rows: list[dict[str, Any]], target: str) -> list[list[str]]:
    satisfying = [row["true_atoms"] for row in rows if row[target]]
    if not satisfying:
        return []
    min_size = min(len(items) for items in satisfying)
    return [items for items in satisfying if len(items) == min_size]


def build_payload() -> dict[str, Any]:
    loaded = {name: load_json(path) for name, path in SOURCE_REPORTS.items()}
    frontier = loaded["theorem_frontier_cut"]
    role_lattice = loaded["role_transfer_minimal_obligation_lattice"]
    component_gap = loaded["component_gap_matrix"]
    anchor = loaded["anchor_h1_classification"]

    atoms = list(frontier["open_atoms"].keys())
    assignments = all_assignments(atoms)
    truth_rows = []
    for idx, assignment in enumerate(assignments):
        values = target_values(assignment)
        true_atoms = true_atom_set(assignment, atoms)
        truth_rows.append(
            {
                "assignment_index": idx,
                "true_atoms": true_atoms,
                "true_atom_count": len(true_atoms),
                **values,
            }
        )

    target_counts = {
        target: sum(1 for row in truth_rows if row[target])
        for target in [
            "bridge_theorem_level_closure",
            "role_transfer_theorem_level_closure",
            "selector_qw2191_closure",
            "toe_closure",
        ]
    }
    minimal_sets = {target: minimal_true_sets(truth_rows, target) for target in target_counts}
    current_assignment = {atom: False for atom in atoms}
    current_targets = target_values(current_assignment)

    summary = {
        "open_atom_count": len(atoms),
        "truth_assignment_count": len(truth_rows),
        "current_assignment_all_false": not any(current_assignment.values()),
        "current_targets_all_false": not any(current_targets.values()),
        "bridge_satisfying_assignment_count": target_counts["bridge_theorem_level_closure"],
        "role_satisfying_assignment_count": target_counts["role_transfer_theorem_level_closure"],
        "selector_satisfying_assignment_count": target_counts["selector_qw2191_closure"],
        "toe_satisfying_assignment_count": target_counts["toe_closure"],
        "bridge_minimal_set_size": len(minimal_sets["bridge_theorem_level_closure"][0]),
        "role_minimal_set_size": len(minimal_sets["role_transfer_theorem_level_closure"][0]),
        "selector_minimal_set_size": len(minimal_sets["selector_qw2191_closure"][0]),
        "toe_minimal_set_size": len(minimal_sets["toe_closure"][0]),
        "toe_minimal_set_equals_frontier_cut": minimal_sets["toe_closure"] == [atoms],
        "role_lattice_minimal_set_inherited": role_lattice["role_transfer_minimal_obligation_summary"]["global_minimal_obligation_size"] == 4,
        "component_gap_sources_still_missing": component_gap["completion_gap_summary"]["strict_dynamic_sources_missing"],
        "anchor_selector_source_still_open": anchor["classification_summary"]["selector_source_remains_open"],
        "bridge_theorem_exported": False,
        "role_transfer_theorem_exported": False,
        "selector_closure_exported": False,
        "toe_closure_claimed": False,
    }

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_THEOREM_FRONTIER_TRUTH_TABLE_CERTIFICATE__NO_CLOSURE_THEOREM",
        "status": "all-2pow7-frontier-assignments-enumerated-current-assignment-closes-no-target",
        "source_reports": {name: rel(path) for name, path in SOURCE_REPORTS.items()},
        "grep_disambiguation": {
            "searched_terms": [
                "truth table",
                "closure assignment",
                "2^7",
                "satisfying assignments",
                "boolean closure",
            ],
            "finding": "No prior strict-completion report enumerated all truth assignments of the seven theorem-frontier atoms; this report adds that finite Boolean audit.",
        },
        "open_atoms": atoms,
        "target_definitions": {
            "bridge_theorem_level_closure": BRIDGE_ATOMS,
            "role_transfer_theorem_level_closure": ROLE_ATOMS,
            "selector_qw2191_closure": SELECTOR_ATOMS,
            "toe_closure": "bridge AND role_transfer AND selector",
        },
        "target_satisfying_assignment_counts": target_counts,
        "target_minimal_true_atom_sets": minimal_sets,
        "current_assignment": current_assignment,
        "current_target_values": current_targets,
        "truth_table_rows": truth_rows,
        "theorem_frontier_truth_table_summary": summary,
        "cross_checks": {
            "source_reports_present": set(loaded) == set(SOURCE_REPORTS),
            "truth_table_size_pass": summary["truth_assignment_count"] == 2 ** summary["open_atom_count"] == 128,
            "current_assignment_no_closure": summary["current_assignment_all_false"] and summary["current_targets_all_false"],
            "minimal_sizes_match_definitions": summary["bridge_minimal_set_size"] == 3 and summary["role_minimal_set_size"] == 4 and summary["selector_minimal_set_size"] == 1 and summary["toe_minimal_set_size"] == 7,
            "toe_minimal_set_matches_frontier_cut": summary["toe_minimal_set_equals_frontier_cut"],
            "limits_preserved": summary["component_gap_sources_still_missing"] and summary["anchor_selector_source_still_open"] and not summary["bridge_theorem_exported"] and not summary["role_transfer_theorem_exported"] and not summary["selector_closure_exported"] and not summary["toe_closure_claimed"],
        },
        "proof_certificate": {
            "grep_step": "rg was used to avoid duplicating the theorem-frontier cut; this report adds the exhaustive Boolean truth table over its seven open atoms.",
            "enumeration_step": "All 2^7=128 truth assignments of the open theorem atoms are enumerated exactly.",
            "target_step": "Bridge closure requires three strict-source atoms; role-transfer closure requires four role/selector atoms; selector closure requires chi11; ToE requires all target predicates.",
            "current_step": "The current assignment has all seven atoms false, so bridge, role-transfer, selector, and ToE targets are all false.",
            "minimal_step": "The minimal true-atom sizes are bridge=3, role=4, selector=1, and ToE=7, matching the frontier-cut report.",
            "scope_step": "This is a Boolean readiness certificate only; it exports no theorem atom, no QW-2191 discharge, and no ToE closure.",
        },
        "hard_limits": [
            "No truth-table assignment is promoted to the current theory state.",
            "No missing theorem atom is exported.",
            "No bridge theorem, role-transfer theorem, selector closure, or ToE closure is claimed.",
            "No QW-2191 selector discharge is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# Strict-completion theorem frontier truth-table certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "All 128 assignments of the seven open theorem atoms are enumerated.",
        "The current all-false assignment closes no target, and ToE closure requires",
        "all seven atoms.",
        "",
        "## Summary",
        "",
    ]
    for key, value in payload["theorem_frontier_truth_table_summary"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Cross-checks", ""])
    for key, value in payload["cross_checks"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Proof certificate", ""])
    for key, value in payload["proof_certificate"].items():
        lines.append(f"- `{key}`: {value}")
    lines.extend(["", "## Hard limits", ""])
    for limit in payload["hard_limits"]:
        lines.append(f"- {limit}")
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    payload = build_payload()
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_markdown(payload)
    print(json.dumps(payload["theorem_frontier_truth_table_summary"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
