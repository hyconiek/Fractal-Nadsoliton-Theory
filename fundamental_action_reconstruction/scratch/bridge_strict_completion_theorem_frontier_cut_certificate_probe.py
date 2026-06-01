#!/usr/bin/env python3
"""Scratch probe: theorem-frontier cut certificate for strict completion.

The bridge ledger now contains finite kernel-comparison witnesses and a role
obligation lattice.  This probe asks the next finite planning question: which
open theorem atoms still form the immediate cut between the current certificate
chain and theorem-level closure?

It builds a small dependency DAG, checks acyclicity by topological sort, computes
transitive prerequisites for the main targets, and reports the minimal open leaf
cut.  It does not export any of those open atoms.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_theorem_frontier_cut_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_theorem_frontier_cut_certificate_report.md"

SOURCE_REPORTS = {
    "finite_bridge_assembly": HERE / "bridge_strict_completion_legacy_to_strict_finite_bridge_assembly_certificate_report.json",
    "symbolic_cancellation": HERE / "bridge_strict_completion_legacy_to_strict_symbolic_cancellation_certificate_report.json",
    "role_transfer_minimal_obligation_lattice": HERE / "bridge_strict_completion_legacy_role_transfer_minimal_obligation_lattice_certificate_report.json",
    "anchor_h1_classification": HERE / "bridge_strict_completion_anchor_h1_generator_classification_certificate_report.json",
    "component_gap_matrix": HERE / "bridge_strict_completion_legacy_to_strict_completion_component_gap_certificate_report.json",
    "closure_plan_dependency": HERE / "bridge_strict_completion_closure_plan_dependency_certificate_report.json",
}

OPEN_ATOMS = [
    "strict_dynamical_source_for_A_P_D",
    "strict_phase_frequency_source",
    "strict_damping_beta_eta_source",
    "chi11_selector_source",
    "alpha_geo_electroweak_role_theorem",
    "beta_tors_strict_role_theorem",
    "beta_power_hierarchy_successor_theorem",
]

CLOSED_ATOMS = [
    "finite_bridge_assembly_witness",
    "symbolic_cancellation_witness",
    "role_transfer_obligation_lattice",
    "anchor_h1_classification",
]

DEPENDENCIES = {
    "finite_bridge_assembly_witness": [],
    "symbolic_cancellation_witness": ["finite_bridge_assembly_witness"],
    "role_transfer_obligation_lattice": ["symbolic_cancellation_witness"],
    "anchor_h1_classification": [],
    "strict_dynamical_source_for_A_P_D": ["symbolic_cancellation_witness"],
    "strict_phase_frequency_source": ["anchor_h1_classification"],
    "strict_damping_beta_eta_source": ["symbolic_cancellation_witness"],
    "chi11_selector_source": ["anchor_h1_classification"],
    "alpha_geo_electroweak_role_theorem": ["role_transfer_obligation_lattice"],
    "beta_tors_strict_role_theorem": ["role_transfer_obligation_lattice"],
    "beta_power_hierarchy_successor_theorem": ["role_transfer_obligation_lattice"],
    "bridge_theorem_level_closure": [
        "strict_dynamical_source_for_A_P_D",
        "strict_phase_frequency_source",
        "strict_damping_beta_eta_source",
    ],
    "role_transfer_theorem_level_closure": [
        "alpha_geo_electroweak_role_theorem",
        "beta_tors_strict_role_theorem",
        "beta_power_hierarchy_successor_theorem",
        "chi11_selector_source",
    ],
    "selector_qw2191_closure": ["chi11_selector_source"],
    "toe_closure": [
        "bridge_theorem_level_closure",
        "role_transfer_theorem_level_closure",
        "selector_qw2191_closure",
    ],
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


def topological_order(graph: dict[str, list[str]]) -> list[str]:
    order: list[str] = []
    temporary: set[str] = set()
    permanent: set[str] = set()

    def visit(node: str) -> None:
        if node in permanent:
            return
        if node in temporary:
            raise ValueError(f"cycle detected at {node}")
        temporary.add(node)
        for dep in graph[node]:
            visit(dep)
        temporary.remove(node)
        permanent.add(node)
        order.append(node)

    for node in graph:
        visit(node)
    return order


def transitive_prerequisites(node: str, graph: dict[str, list[str]]) -> set[str]:
    result: set[str] = set()
    for dep in graph[node]:
        result.add(dep)
        result |= transitive_prerequisites(dep, graph)
    return result


def build_payload() -> dict[str, Any]:
    loaded = {name: load_json(path) for name, path in SOURCE_REPORTS.items()}
    finite_bridge = loaded["finite_bridge_assembly"]
    symbolic = loaded["symbolic_cancellation"]
    role_lattice = loaded["role_transfer_minimal_obligation_lattice"]
    anchor = loaded["anchor_h1_classification"]
    component_gap = loaded["component_gap_matrix"]
    closure_plan = loaded["closure_plan_dependency"]

    closed_status = {
        "finite_bridge_assembly_witness": finite_bridge["finite_bridge_assembly_summary"]["assembled_map_reconstructs_strict_kernel"],
        "symbolic_cancellation_witness": symbolic["symbolic_cancellation_summary"]["symbolic_cancellation_formula_exported"],
        "role_transfer_obligation_lattice": role_lattice["role_transfer_minimal_obligation_summary"]["global_minimal_obligation_count"] == 1,
        "anchor_h1_classification": anchor["classification_summary"]["left_anchor_is_c0_gauge_fix_not_c1_generator"],
    }
    open_status = {atom: False for atom in OPEN_ATOMS}
    topo = topological_order(DEPENDENCIES)

    target_rows = []
    all_open_leaf_cut: set[str] = set()
    for target in TARGETS:
        prereqs = transitive_prerequisites(target, DEPENDENCIES)
        open_prereqs = sorted(atom for atom in prereqs if atom in OPEN_ATOMS)
        closed_prereqs = sorted(atom for atom in prereqs if atom in CLOSED_ATOMS)
        all_open_leaf_cut.update(open_prereqs)
        target_rows.append(
            {
                "target": target,
                "direct_dependencies": DEPENDENCIES[target],
                "closed_prerequisites": closed_prereqs,
                "open_leaf_cut": open_prereqs,
                "open_leaf_cut_size": len(open_prereqs),
                "target_exported_now": False,
            }
        )

    summary = {
        "dag_node_count": len(DEPENDENCIES),
        "topological_order": topo,
        "dag_is_acyclic": len(topo) == len(DEPENDENCIES),
        "closed_atom_count": sum(closed_status.values()),
        "open_atom_count": len(OPEN_ATOMS),
        "all_closed_atoms_certified": all(closed_status.values()),
        "all_open_atoms_still_missing": not any(open_status.values()),
        "minimal_open_leaf_cut_for_toe": next(row for row in target_rows if row["target"] == "toe_closure")["open_leaf_cut"],
        "minimal_open_leaf_cut_for_toe_size": next(row for row in target_rows if row["target"] == "toe_closure")["open_leaf_cut_size"],
        "union_open_leaf_cut": sorted(all_open_leaf_cut),
        "component_gap_sources_still_missing": component_gap["completion_gap_summary"]["strict_dynamic_sources_missing"],
        "closure_plan_keeps_toe_open": closure_plan["closure_plan_summary"]["toe_closure_not_claimed"],
        "bridge_theorem_exported": False,
        "role_transfer_theorem_exported": False,
        "selector_closure_exported": False,
        "toe_closure_claimed": False,
    }

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_THEOREM_FRONTIER_CUT_CERTIFICATE__NO_CLOSURE_THEOREM",
        "status": "theorem-frontier-dag-acyclic-open-leaf-cut-enumerated-no-closure-exported",
        "source_reports": {name: rel(path) for name, path in SOURCE_REPORTS.items()},
        "grep_disambiguation": {
            "searched_terms": [
                "theorem readiness",
                "bridge closure readiness",
                "frontier cut",
                "minimal cut",
                "source obligation",
                "selector obligation",
            ],
            "finding": "Existing reports certify finite bridge and role-obligation objects; this report computes the remaining theorem-frontier cut instead of redoing those certificates.",
        },
        "closed_atoms": closed_status,
        "open_atoms": open_status,
        "dependency_graph": DEPENDENCIES,
        "target_frontier_rows": target_rows,
        "theorem_frontier_cut_summary": summary,
        "cross_checks": {
            "source_reports_present": set(loaded) == set(SOURCE_REPORTS),
            "dag_acyclic_and_closed_atoms_loaded": summary["dag_is_acyclic"] and summary["all_closed_atoms_certified"],
            "open_atoms_missing": summary["all_open_atoms_still_missing"] and summary["open_atom_count"] == len(OPEN_ATOMS),
            "toe_cut_contains_all_source_and_role_atoms": set(summary["minimal_open_leaf_cut_for_toe"]) == set(OPEN_ATOMS),
            "prior_closure_limits_preserved": summary["component_gap_sources_still_missing"] and summary["closure_plan_keeps_toe_open"] and not summary["bridge_theorem_exported"] and not summary["role_transfer_theorem_exported"] and not summary["selector_closure_exported"] and not summary["toe_closure_claimed"],
        },
        "proof_certificate": {
            "grep_step": "rg was used to distinguish this theorem-frontier cut from existing closure-plan and role-obligation reports.",
            "dag_step": f"The dependency graph has {len(DEPENDENCIES)} nodes and is acyclic by topological sort.",
            "closed_step": "Finite bridge assembly, symbolic cancellation, role obligation lattice, and anchor/H1 classification are treated as closed certificate atoms.",
            "open_cut_step": "The ToE frontier cut consists exactly of seven missing theorem atoms: three strict source atoms, three legacy role atoms, and chi11 selector/source.",
            "scope_step": "This is a readiness/cut certificate only; it exports no bridge theorem, role-transfer theorem, selector closure, QW-2191 discharge, or ToE closure.",
        },
        "hard_limits": [
            "No missing theorem atom is exported by this frontier-cut certificate.",
            "No bridge theorem-level closure is claimed.",
            "No legacy role-transfer theorem is claimed.",
            "No selector/QW-2191 closure is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# Strict-completion theorem frontier cut certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "A finite dependency DAG separates already certified finite witnesses from",
        "open theorem atoms.  The ToE target remains blocked by the seven-atom open",
        "leaf cut; none of those atoms is exported here.",
        "",
        "## Summary",
        "",
    ]
    for key, value in payload["theorem_frontier_cut_summary"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Target rows", ""])
    for row in payload["target_frontier_rows"]:
        lines.append(f"- `{row['target']}` open cut: `{row['open_leaf_cut']}`")
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
    print(json.dumps(payload["theorem_frontier_cut_summary"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
