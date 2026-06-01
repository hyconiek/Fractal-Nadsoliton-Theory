#!/usr/bin/env python3
"""Scratch probe: theorem-frontier atom influence certificate.

The truth-table certificate enumerates all 2^7 open-frontier assignments.  This
probe takes the next finite step: for each open theorem atom and each closure
target, compute the Boolean swing/critical count, i.e. the number of assignments
where toggling that atom from false to true changes the target from false to
true.

The result is a prioritization/audit metric only.  It does not export any
missing theorem atom and does not close the bridge, role-transfer, selector, or
ToE targets.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_theorem_frontier_atom_influence_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_theorem_frontier_atom_influence_certificate_report.md"

SOURCE_REPORTS = {
    "theorem_frontier_truth_table": HERE / "bridge_strict_completion_theorem_frontier_truth_table_certificate_report.json",
    "theorem_frontier_cut": HERE / "bridge_strict_completion_theorem_frontier_cut_certificate_report.json",
    "role_transfer_minimal_obligation_lattice": HERE / "bridge_strict_completion_legacy_role_transfer_minimal_obligation_lattice_certificate_report.json",
}

TARGETS = [
    "bridge_theorem_level_closure",
    "role_transfer_theorem_level_closure",
    "selector_qw2191_closure",
    "toe_closure",
]


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def target_values(assignment: dict[str, bool], target_defs: dict[str, Any]) -> dict[str, bool]:
    bridge = all(assignment[atom] for atom in target_defs["bridge_theorem_level_closure"])
    role = all(assignment[atom] for atom in target_defs["role_transfer_theorem_level_closure"])
    selector = all(assignment[atom] for atom in target_defs["selector_qw2191_closure"])
    toe = bridge and role and selector
    return {
        "bridge_theorem_level_closure": bridge,
        "role_transfer_theorem_level_closure": role,
        "selector_qw2191_closure": selector,
        "toe_closure": toe,
    }


def assignments_with_atom_false(atoms: list[str], atom: str) -> list[dict[str, bool]]:
    other_atoms = [item for item in atoms if item != atom]
    rows = []
    for bits in itertools.product([False, True], repeat=len(other_atoms)):
        assignment = dict(zip(other_atoms, bits))
        assignment[atom] = False
        rows.append({name: assignment[name] for name in atoms})
    return rows


def critical_count(atoms: list[str], atom: str, target: str, target_defs: dict[str, Any]) -> int:
    count = 0
    for before in assignments_with_atom_false(atoms, atom):
        after = dict(before)
        after[atom] = True
        if not target_values(before, target_defs)[target] and target_values(after, target_defs)[target]:
            count += 1
    return count


def target_memberships(atom: str, target_defs: dict[str, Any]) -> list[str]:
    memberships = []
    for target in ["bridge_theorem_level_closure", "role_transfer_theorem_level_closure", "selector_qw2191_closure"]:
        if atom in target_defs[target]:
            memberships.append(target)
    memberships.append("toe_closure")
    return memberships


def build_payload() -> dict[str, Any]:
    loaded = {name: load_json(path) for name, path in SOURCE_REPORTS.items()}
    truth_table = loaded["theorem_frontier_truth_table"]
    frontier_cut = loaded["theorem_frontier_cut"]
    role_lattice = loaded["role_transfer_minimal_obligation_lattice"]

    atoms = truth_table["open_atoms"]
    target_defs = truth_table["target_definitions"]
    rows = []
    for atom in atoms:
        counts = {target: critical_count(atoms, atom, target, target_defs) for target in TARGETS}
        total = sum(counts.values())
        rows.append(
            {
                "atom": atom,
                "target_memberships": target_memberships(atom, target_defs),
                "bridge_critical_count": counts["bridge_theorem_level_closure"],
                "role_transfer_critical_count": counts["role_transfer_theorem_level_closure"],
                "selector_qw2191_critical_count": counts["selector_qw2191_closure"],
                "toe_critical_count": counts["toe_closure"],
                "total_critical_count": total,
                "normalized_total_critical_count": total / (len(TARGETS) * (2 ** (len(atoms) - 1))),
                "exported_now": False,
            }
        )

    rows_by_atom = {row["atom"]: row for row in rows}
    ordered = sorted(rows, key=lambda item: (-item["total_critical_count"], item["atom"]))
    top_count = ordered[0]["total_critical_count"]
    top_atoms = [row["atom"] for row in ordered if row["total_critical_count"] == top_count]

    summary = {
        "open_atom_count": len(atoms),
        "target_count": len(TARGETS),
        "critical_pair_universe_per_atom_target": 2 ** (len(atoms) - 1),
        "total_critical_count_all_atoms": sum(row["total_critical_count"] for row in rows),
        "top_influence_atoms": top_atoms,
        "top_influence_total_critical_count": top_count,
        "chi11_selector_source_is_unique_top_influence": top_atoms == ["chi11_selector_source"],
        "chi11_selector_source_total_critical_count": rows_by_atom["chi11_selector_source"]["total_critical_count"],
        "bridge_source_atoms_tie_at_17": all(rows_by_atom[atom]["total_critical_count"] == 17 for atom in target_defs["bridge_theorem_level_closure"]),
        "role_only_atoms_tie_at_9": all(rows_by_atom[atom]["total_critical_count"] == 9 for atom in ["alpha_geo_electroweak_role_theorem", "beta_tors_strict_role_theorem", "beta_power_hierarchy_successor_theorem"]),
        "each_atom_is_toe_critical_once": all(row["toe_critical_count"] == 1 for row in rows),
        "truth_table_current_assignment_closes_no_target": truth_table["theorem_frontier_truth_table_summary"]["current_targets_all_false"],
        "frontier_cut_open_atoms_match": atoms == list(frontier_cut["open_atoms"].keys()),
        "role_lattice_atoms_remain_missing": role_lattice["role_transfer_minimal_obligation_summary"]["all_atoms_missing_now"],
        "bridge_theorem_exported": False,
        "role_transfer_theorem_exported": False,
        "selector_closure_exported": False,
        "toe_closure_claimed": False,
    }

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_THEOREM_FRONTIER_ATOM_INFLUENCE_CERTIFICATE__NO_THEOREM_EXPORT",
        "status": "boolean-swing-criticality-computed-for-seven-open-frontier-atoms-no-closure",
        "source_reports": {name: rel(path) for name, path in SOURCE_REPORTS.items()},
        "grep_disambiguation": {
            "searched_terms": [
                "Banzhaf",
                "critical count",
                "atom influence",
                "Boolean derivative",
                "bottleneck score",
            ],
            "finding": "No prior strict-completion report computed Boolean swing/criticality counts for each open theorem-frontier atom; this report adds that finite influence audit.",
        },
        "open_atoms": atoms,
        "targets": TARGETS,
        "atom_influence_rows": rows,
        "atom_influence_order": [row["atom"] for row in ordered],
        "theorem_frontier_atom_influence_summary": summary,
        "cross_checks": {
            "source_reports_present": set(loaded) == set(SOURCE_REPORTS),
            "atom_set_matches_truth_table_and_cut": summary["frontier_cut_open_atoms_match"] and atoms == truth_table["open_atoms"],
            "critical_counts_match_closed_forms": rows_by_atom["chi11_selector_source"]["total_critical_count"] == 73 and all(rows_by_atom[atom]["total_critical_count"] == 17 for atom in target_defs["bridge_theorem_level_closure"]) and summary["role_only_atoms_tie_at_9"],
            "top_influence_is_chi11_selector_source": summary["chi11_selector_source_is_unique_top_influence"],
            "toe_criticality_is_uniform": summary["each_atom_is_toe_critical_once"],
            "limits_preserved": summary["truth_table_current_assignment_closes_no_target"] and summary["role_lattice_atoms_remain_missing"] and not summary["bridge_theorem_exported"] and not summary["role_transfer_theorem_exported"] and not summary["selector_closure_exported"] and not summary["toe_closure_claimed"],
        },
        "proof_certificate": {
            "grep_step": "rg was used to avoid duplicating an existing influence/power-index audit; none existed for the strict-completion theorem frontier.",
            "definition_step": "For an atom a and target T, the critical count enumerates assignments with a=false where flipping a to true changes T from false to true.",
            "enumeration_step": "Each atom-target pair checks 2^6=64 assignments, so the audit covers 7*4*64 finite swing cases.",
            "influence_step": "chi11_selector_source has total critical count 73 because it is critical for selector in 64 assignments, role transfer in 8 assignments, and ToE in 1 assignment.",
            "bridge_step": "Each strict bridge source atom has total critical count 17: 16 bridge swings plus the single ToE all-other-true swing.",
            "role_step": "The three role-only atoms have total critical count 9: 8 role-transfer swings plus the single ToE all-other-true swing.",
            "scope_step": "This ranks open theorem atoms by Boolean criticality only; it exports no atom and proves no bridge, role-transfer, selector, or ToE closure.",
        },
        "hard_limits": [
            "No influence score is promoted to a theorem source.",
            "No missing theorem atom is exported.",
            "No bridge theorem, role-transfer theorem, selector closure, or ToE closure is claimed.",
            "No QW-2191 selector discharge is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# Strict-completion theorem frontier atom-influence certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "Boolean swing/criticality counts are computed for every open theorem atom.",
        "The unique top logical bottleneck is `chi11_selector_source`, but this is",
        "only a finite influence audit and exports no selector source theorem.",
        "",
        "## Summary",
        "",
    ]
    for key, value in payload["theorem_frontier_atom_influence_summary"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Atom influence rows", ""])
    for row in payload["atom_influence_rows"]:
        lines.append(
            f"- `{row['atom']}`: total=`{row['total_critical_count']}`, "
            f"bridge=`{row['bridge_critical_count']}`, role=`{row['role_transfer_critical_count']}`, "
            f"selector=`{row['selector_qw2191_critical_count']}`, ToE=`{row['toe_critical_count']}`"
        )
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
    print(json.dumps(payload["theorem_frontier_atom_influence_summary"], indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
