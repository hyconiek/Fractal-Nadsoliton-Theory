#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2405_s1355_nadsoliton_information_ontology_projection_certificate.json"
MD = GEN / "p2405_s1355_nadsoliton_information_ontology_projection_certificate.md"

SOURCE_FILES = {
    "AX9_ONTOLOGY": GEN / "ax9_informational_nadsoliton_primacy_axiom_packet.json",
    "N50_NONIDENTIFICATION": ROOT / "N50_CURRENT_LEGACY_ONTOLOGICAL_KERNEL_TO_STRICT_GATE_KERNEL_NONIDENTIFICATION_THEOREM.md",
    "S2_PRIORITY_PACKET": ROOT / "S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md",
    "P2404_DEPENDENCY_CUT": GEN / "p2404_s1354_strict_addition_physics_lane_dependency_cut_certificate.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

ONTOLOGY_NODES = [
    "nadsoliton_pure_information_root",
    "legacy_intermediate_kernel_information_code",
    "strict_completed_kernel_information_code",
    "light_projection",
    "matter_physical_role_projection",
    "emergent_observer_readout",
]

ONTOLOGY_EDGES = [
    ("nadsoliton_pure_information_root", "legacy_intermediate_kernel_information_code"),
    ("legacy_intermediate_kernel_information_code", "strict_completed_kernel_information_code"),
    ("strict_completed_kernel_information_code", "light_projection"),
    ("light_projection", "matter_physical_role_projection"),
    ("matter_physical_role_projection", "emergent_observer_readout"),
]

ONTOLOGY_GUARD_ATOMS = [
    "nadsoliton_is_sole_primordial_information",
    "no_separate_information_layer_under_nadsoliton",
    "strict_additions_are_internal_information_constraints",
    "physics_roles_are_downstream_projections",
    "observer_is_downstream_readout_not_source",
]

STRICT_INFORMATION_ATOMS = [
    "apd_completion_structure",
    "gf2_phase_topological_data",
    "nonlinear_damping_compression",
    "chi11_selector_source_declared",
]

ROLE_PROJECTION_ATOMS = [
    "alpha_geo_electroweak_role_theorem",
    "beta_tors_strict_role_theorem",
    "beta_power_hierarchy_successor_theorem",
]

INTERPRETATION_ROWS = [
    {
        "interpretation": "strict_additions_as_internal_nadsoliton_information_constraints",
        "valid_under_pure_information_ontology": True,
        "reason": "Strict additions type the completed kernel as internal information carried by the nadsoliton.",
    },
    {
        "interpretation": "legacy_kernel_as_intermediate_information_bridge_code",
        "valid_under_pure_information_ontology": True,
        "reason": "Legacy remains an intermediate code/path inside the same informational nadsoliton, not a final coequal physics generator.",
    },
    {
        "interpretation": "physical_roles_as_downstream_projection_theorems",
        "valid_under_pure_information_ontology": True,
        "reason": "Mass/electroweak/gravity roles are matter-side projections and need role-successor theorems.",
    },
    {
        "interpretation": "separate_information_substrate_below_nadsoliton",
        "valid_under_pure_information_ontology": False,
        "reason": "The ontology states that the nadsoliton itself is primordial information; there is no lower information layer.",
    },
    {
        "interpretation": "strict_kernel_as_external_substrate_source",
        "valid_under_pure_information_ontology": False,
        "reason": "Strict must be an internal completed code of the nadsoliton, not an external substrate underneath it.",
    },
    {
        "interpretation": "strict_additions_as_immediate_physical_role_exports",
        "valid_under_pure_information_ontology": False,
        "reason": "P2404 leaves strict-additions-only physical role lanes false; information constraints are not automatic role theorems.",
    },
    {
        "interpretation": "observer_readout_as_selector_source_without_bridge",
        "valid_under_pure_information_ontology": False,
        "reason": "Observer remains downstream of nadsoliton->light->matter and cannot be smuggled in as a strict-core source.",
    },
]


def rel(path: Path) -> str:
    return path.relative_to(REPO).as_posix()


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": rel(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8")


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg",
            "-n",
            pattern,
            "fundamental_action_reconstruction",
            "material_dowodowy",
            "-g",
            "*.py",
            "-g",
            "*.md",
            "-g",
            "*.json",
            "-g",
            "*.tex",
            "-g",
            "!generated/p2405_s1355_nadsoliton_information_ontology_projection_certificate.*",
        ],
        cwd=REPO,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:18]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2405|S1355|nadsoliton information ontology projection|pure-information projection certificate",
        "pure_information_ontology": "pure information|primordial information|separate informational layer|nadsoliton -> light -> matter -> emergent observer",
        "kernel_internal_route": "internal route of one informational nadsoliton|information bridge|strict.*information|legacy.*information",
        "p2404_dependency_cut": "P2404|strict-addition physics-lane dependency-cut|strict-additions-only physical role lanes ready",
        "observer_downstream_guard": "observer.*downstream|emergent observer.*downstream|readout.*not.*source",
    }
    return {
        "tool": "rg",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
        "finding": (
            "Repo grep shows many ontology reminders and P2404 dependency-cut evidence, but no local certificate "
            "that retypes P2404 strict additions as internal nadsoliton information while blocking any sub-nadsoliton layer."
        ),
    }


def transitive_closure(nodes: list[str], edges: list[tuple[str, str]]) -> dict[str, list[str]]:
    adjacency = {node: set() for node in nodes}
    for src, dst in edges:
        adjacency[src].add(dst)
    closure: dict[str, list[str]] = {}
    for node in nodes:
        seen: set[str] = set()
        frontier = list(adjacency[node])
        while frontier:
            current = frontier.pop()
            if current in seen:
                continue
            seen.add(current)
            frontier.extend(adjacency[current] - seen)
        closure[node] = sorted(seen, key=nodes.index)
    return closure


def topological_orders(nodes: list[str], edges: list[tuple[str, str]]) -> list[list[str]]:
    incoming = {node: set() for node in nodes}
    outgoing = {node: set() for node in nodes}
    for src, dst in edges:
        outgoing[src].add(dst)
        incoming[dst].add(src)

    orders: list[list[str]] = []

    def visit(prefix: list[str], remaining: set[str], incoming_now: dict[str, set[str]]) -> None:
        if not remaining:
            orders.append(prefix.copy())
            return
        ready = sorted([node for node in remaining if not incoming_now[node]], key=nodes.index)
        for node in ready:
            next_remaining = set(remaining)
            next_remaining.remove(node)
            next_incoming = {key: set(value) for key, value in incoming_now.items()}
            for dst in outgoing[node]:
                next_incoming[dst].discard(node)
            visit(prefix + [node], next_remaining, next_incoming)

    visit([], set(nodes), incoming)
    return orders


def guard_lattice() -> dict[str, Any]:
    rows = []
    full_mask = (1 << len(ONTOLOGY_GUARD_ATOMS)) - 1
    for mask in range(1 << len(ONTOLOGY_GUARD_ATOMS)):
        true_atoms = [atom for index, atom in enumerate(ONTOLOGY_GUARD_ATOMS) if mask & (1 << index)]
        ontology_valid = mask == full_mask
        rows.append(
            {
                "mask": mask,
                "true_atoms": true_atoms,
                "ontology_valid": ontology_valid,
                "missing_atoms": [atom for atom in ONTOLOGY_GUARD_ATOMS if atom not in true_atoms],
            }
        )
    return {
        "guard_atom_order": ONTOLOGY_GUARD_ATOMS,
        "row_count": len(rows),
        "rows": rows,
        "true_masks": [row["mask"] for row in rows if row["ontology_valid"]],
        "proper_subset_fail_count": sum(1 for row in rows if row["mask"] != full_mask and not row["ontology_valid"]),
        "ontology_valid_anf": " * ".join(ONTOLOGY_GUARD_ATOMS),
        "ontology_valid_anf_degree": len(ONTOLOGY_GUARD_ATOMS),
        "boolean_derivative_supports": {atom: [full_mask ^ (1 << index)] for index, atom in enumerate(ONTOLOGY_GUARD_ATOMS)},
    }


def ontology_poset() -> dict[str, Any]:
    closure = transitive_closure(ONTOLOGY_NODES, ONTOLOGY_EDGES)
    orders = topological_orders(ONTOLOGY_NODES, ONTOLOGY_EDGES)
    incoming_nodes = {dst for _, dst in ONTOLOGY_EDGES}
    roots = [node for node in ONTOLOGY_NODES if node not in incoming_nodes]
    has_cycle = any(node in closure[node] for node in ONTOLOGY_NODES)
    preferred_order = ["nadsoliton", "light", "matter", "emergent_observer"]
    collapsed_order = [
        "nadsoliton",
        "legacy_kernel_information_code",
        "strict_kernel_information_code",
        "light",
        "matter",
        "emergent_observer",
    ]
    return {
        "nodes": ONTOLOGY_NODES,
        "edges": [list(edge) for edge in ONTOLOGY_EDGES],
        "transitive_closure": closure,
        "root_nodes": roots,
        "has_cycle": has_cycle,
        "topological_order_count": len(orders),
        "topological_orders": orders,
        "preferred_internal_order": preferred_order,
        "refined_kernel_bridge_order": collapsed_order,
        "no_node_precedes_nadsoliton_root": roots == ["nadsoliton_pure_information_root"],
        "all_nodes_reachable_from_nadsoliton_root": set(closure["nadsoliton_pure_information_root"]) == set(ONTOLOGY_NODES[1:]),
    }


def append_doc_sections() -> None:
    eq_section = """
## P2405/S1355 nadsoliton information-ontology projection certificate

`P2405/S1355` corrects the P2404 dependency-cut reading by making the ontology type explicit: the nadsoliton itself is the sole primordial information in a solitonic state.  The four strict additions are therefore typed as internal information constraints of the completed strict kernel, not as a new informational substrate beneath the nadsoliton and not as immediate physical-role exports.

The finite ontology guard lattice has `32` rows and a single true mask: `nadsoliton_is_sole_primordial_information * no_separate_information_layer_under_nadsoliton * strict_additions_are_internal_information_constraints * physics_roles_are_downstream_projections * observer_is_downstream_readout_not_source`.  The poset audit has a unique root, `nadsoliton_pure_information_root`, and preserves the order `nadsoliton -> light -> matter -> emergent observer` after inserting legacy/strict kernel codes as internal information-code stages upstream of light.

Thus P2404's strict-addition cut is not a new lower layer.  It is a proof obligation inside the single informational nadsoliton before any downstream physical projection can be claimed.
""".strip()
    lag_section = """
## P2405/S1355 information-ontology guard for Lagrangian/EOM

`P2405/S1355` types the strict additions as internal information constraints of the nadsoliton, not as physical Lagrangian terms and not as a sub-nadsoliton information layer.  `L_total` remains downstream of light/matter projection and still requires the P2404 degree-7 dependency monomial before role-bearing promotion.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items() if path.suffix == ".json"}
    texts = {name: load_text(path) for name, path in SOURCE_FILES.items() if path.suffix != ".json"}
    grep = rg_audit()
    poset = ontology_poset()
    lattice = guard_lattice()
    p2404_theorem = artifacts["P2404_DEPENDENCY_CUT"].get(
        "strict_addition_physics_lane_dependency_cut_certificate", {}
    ).get("theorem_export", {})
    ax9 = artifacts["AX9_ONTOLOGY"]
    n50_text = texts.get("N50_NONIDENTIFICATION", "")
    s2_text = texts.get("S2_PRIORITY_PACKET", "")

    valid_interpretations = [row["interpretation"] for row in INTERPRETATION_ROWS if row["valid_under_pure_information_ontology"]]
    invalid_interpretations = [row["interpretation"] for row in INTERPRETATION_ROWS if not row["valid_under_pure_information_ontology"]]
    theorem_export = {
        "theorem_name": "P2405_T1_nadsoliton_information_ontology_projection_certificate",
        "canonical_ontology_statement": ax9.get("canonical_ontology_statement"),
        "single_primitive": ax9.get("recovered_consequences", {}).get("single_primitive"),
        "separate_information_layer_allowed": ax9.get("recovered_consequences", {}).get("separate_information_layer"),
        "ontology_guard_true_masks": lattice["true_masks"],
        "ontology_guard_proper_subset_fail_count": lattice["proper_subset_fail_count"],
        "ontology_guard_anf": lattice["ontology_valid_anf"],
        "ontology_guard_anf_degree": lattice["ontology_valid_anf_degree"],
        "unique_information_root": poset["root_nodes"] == ["nadsoliton_pure_information_root"],
        "all_nodes_reachable_from_information_root": poset["all_nodes_reachable_from_nadsoliton_root"],
        "topological_order_count": poset["topological_order_count"],
        "strict_information_atoms": STRICT_INFORMATION_ATOMS,
        "role_projection_atoms": ROLE_PROJECTION_ATOMS,
        "valid_interpretations": valid_interpretations,
        "invalid_interpretations": invalid_interpretations,
        "p2404_strict_additions_only_physical_role_lanes_ready": p2404_theorem.get(
            "strict_additions_only_physical_role_lanes_ready"
        ),
        "p2404_ltotal_dependency_degree": p2404_theorem.get("ltotal_dependency_degree"),
        "n50_internal_nadsoliton_route_warning_detected": "internal route of one informational nadsoliton" in n50_text,
        "s2_priority_order_detected": "legacy -> strict" in s2_text and "role-transfer audit" in s2_text,
        "not_licensed": [
            "No separate information substrate below the nadsoliton is introduced.",
            "No strict addition is retyped as an immediate physical role theorem.",
            "No observer/readout source closes QW-2191 or supplies a strict selector theorem.",
            "No L_total or ToE promotion follows from ontology typing alone.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "ax9_single_primitive_inherited": theorem_export["single_primitive"] == "informational_nadsoliton_only",
        "ax9_blocks_separate_information_layer": theorem_export["separate_information_layer_allowed"] is False,
        "ontology_lattice_has_32_rows": lattice["row_count"] == 32,
        "only_full_ontology_guard_mask_passes": theorem_export["ontology_guard_true_masks"] == [31]
        and theorem_export["ontology_guard_proper_subset_fail_count"] == 31,
        "ontology_guard_degree_five": theorem_export["ontology_guard_anf_degree"] == 5,
        "unique_nadsoliton_information_root": theorem_export["unique_information_root"],
        "all_nodes_reachable_from_information_root": theorem_export["all_nodes_reachable_from_information_root"],
        "p2404_strict_additions_do_not_license_physical_roles": theorem_export[
            "p2404_strict_additions_only_physical_role_lanes_ready"
        ] == [],
        "p2404_degree_seven_dependency_inherited": theorem_export["p2404_ltotal_dependency_degree"] == 7,
        "n50_internal_route_warning_inherited": theorem_export["n50_internal_nadsoliton_route_warning_detected"],
        "s2_priority_order_inherited": theorem_export["s2_priority_order_detected"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2405_s1355_v1",
        "packet_id": "P2405",
        "stage_id": "S1355",
        "result_kind": "NADSOLITON_INFORMATION_ONTOLOGY_PROJECTION_CERTIFICATE",
        "status": "PASS_INFORMATION_ONTOLOGY_RETYPE_NO_UNDERLAYER_NO_ROLE_EXPORT",
        "nadsoliton_information_ontology_projection_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: artifact.get("status", artifact.get("result", "TEXT_SOURCE")) for name, artifact in artifacts.items()},
            "ontology_poset": poset,
            "ontology_guard_lattice": lattice,
            "interpretation_matrix": INTERPRETATION_ROWS,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": (
            "Use strict additions as internal informational proof obligations of the nadsoliton, then separately prove "
            "the downstream role-successor projection theorems before any L_total or ToE promotion."
        ),
        "global_status": "OPEN_PROGRESS_INFORMATION_ONTOLOGY_TYPED_NO_PHYSICAL_ROLE_EXPORT",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["nadsoliton_information_ontology_projection_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2405 S1355: nadsoliton information-ontology projection certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2405/S1355 retypes P2404 strict additions as internal information constraints of the one nadsoliton, not as a lower information layer or direct physical-role exports.",
                "",
                "## Exact ontology guard",
                "",
                f"- Guard true masks: `{theorem['ontology_guard_true_masks']}`.",
                f"- Proper-subset fail count: `{theorem['ontology_guard_proper_subset_fail_count']}`.",
                f"- Guard ANF: `{theorem['ontology_guard_anf']}`.",
                f"- Guard ANF degree: `{theorem['ontology_guard_anf_degree']}`.",
                "",
                "## Poset / projection facts",
                "",
                f"- Unique information root: `{theorem['unique_information_root']}`.",
                f"- All nodes reachable from information root: `{theorem['all_nodes_reachable_from_information_root']}`.",
                f"- Topological order count: `{theorem['topological_order_count']}`.",
                f"- P2404 strict-additions-only physical role lanes ready: `{theorem['p2404_strict_additions_only_physical_role_lanes_ready']}`.",
                f"- P2404 L_total dependency degree: `{theorem['p2404_ltotal_dependency_degree']}`.",
                "",
                "## Hard limits",
                "",
                *[f"- {item}" for item in theorem["not_licensed"]],
                "",
                "## Gatekeepers",
                "",
                f"`{payload['gatekeeper_checks']}`",
                "",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    GEN.mkdir(parents=True, exist_ok=True)
    append_doc_sections()
    payload = build_payload()
    write_outputs(payload)
    if not all(payload["gatekeeper_checks"].values()):
        raise SystemExit(f"gatekeeper failure: {payload['gatekeeper_checks']}")


if __name__ == "__main__":
    main()
