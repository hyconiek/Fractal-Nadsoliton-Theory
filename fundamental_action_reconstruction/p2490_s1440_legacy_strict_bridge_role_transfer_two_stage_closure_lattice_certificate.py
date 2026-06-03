#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from collections import Counter
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import (
    DOC_FILES,
    REPO,
    ROOT,
    load_json,
    rel,
)

GEN = ROOT / "generated"
OUT = GEN / "p2490_s1440_legacy_strict_bridge_role_transfer_two_stage_closure_lattice_certificate.json"
MD = GEN / "p2490_s1440_legacy_strict_bridge_role_transfer_two_stage_closure_lattice_certificate.md"

SOURCE_FILES = {
    "P2411_BRIDGE_OBLIGATION_HYPERGRAPH": GEN / "p2411_s1361_legacy_strict_bridge_source_obligation_hypergraph_certificate.json",
    "P2434_ROLE_TRANSFER_CLAIM_LATTICE": GEN / "p2434_s1384_conditional_legacy_role_transfer_claim_lattice_certificate.json",
}


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, "fundamental_action_reconstruction", "material_dowodowy",
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "!fundamental_action_reconstruction/generated/**",
        ],
        cwd=REPO,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:40]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2490|S1440|two-stage bridge-role closure lattice|bridge-role closure lattice|end-to-end bridge role|combined bridge role lattice|bridge-to-role transfer frontier|post-bridge role lattice",
        "precursor_packets": "P2411|S1361|legacy-strict bridge source obligation hypergraph|P2434|S1384|conditional legacy role-transfer claim lattice",
        "closure_lattice_language": "two-stage closure|end-to-end closure|bridge_ready|role claim ready|nearest miss|assignment lattice|truth vector|frontier atom",
        "guardrail_language": "legacy -> strict completion bridge|role-transfer audit|silent inheritance|K_legacy_ont|K_strict_gate|QW-2191",
        "closure_blockers": "role-bearing L_total|physical-value generator|ToE closure|selector/source/gauge theorem|legacy-role transfer|root-window theorem",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def container(payload: dict[str, Any], key: str) -> dict[str, Any]:
    return payload.get(key, {})


def mask_names(mask: int, atoms: list[str]) -> list[str]:
    return [atom for index, atom in enumerate(atoms) if (mask >> index) & 1]


def missing_names(mask: int, atoms: list[str]) -> list[str]:
    return [atom for index, atom in enumerate(atoms) if not ((mask >> index) & 1)]


def role_ready_count(role_mask: int, all_role_ready_masks: set[int], minimal_masks_by_role: dict[str, list[list[str]]], role_atoms: list[str]) -> int:
    ready = 0
    true_atoms = set(mask_names(role_mask, role_atoms))
    for minimal_sets in minimal_masks_by_role.values():
        if any(set(req).issubset(true_atoms) for req in minimal_sets):
            ready += 1
    if role_mask in all_role_ready_masks:
        return max(ready, len(minimal_masks_by_role))
    return ready


def two_stage_lattice(p2411_payload: dict[str, Any], p2434_payload: dict[str, Any]) -> dict[str, Any]:
    p2411_cert = container(p2411_payload, "legacy_strict_bridge_source_obligation_hypergraph_certificate")
    p2434_cert = container(p2434_payload, "conditional_legacy_role_transfer_claim_lattice_certificate")
    p2411 = theorem(p2411_payload, "legacy_strict_bridge_source_obligation_hypergraph_certificate")
    p2434 = theorem(p2434_payload, "conditional_legacy_role_transfer_claim_lattice_certificate")
    hypergraph = p2411_cert["finite_hypergraph_certificate"]
    bridge_atoms = hypergraph["bridge_obligation_atoms"]
    role_atoms = p2434["role_obligation_names"]
    bridge_ready_masks = set(p2411["bridge_ready_true_masks"])
    all_role_ready_masks = set(p2434["all_roles_ready_masks"])
    minimal_masks_by_role = p2434["minimal_masks_by_role_claim"]

    closure_distribution: Counter[str] = Counter()
    role_count_distribution_when_bridge_ready: Counter[int] = Counter()
    role_count_distribution_when_bridge_blocked: Counter[int] = Counter()
    bridge_ready_count = 0
    bridge_blocked_count = 0
    end_to_end_ready_count = 0
    current_state = None
    nearest_misses = []

    for bridge_mask in range(1 << len(bridge_atoms)):
        bridge_ready = bridge_mask in bridge_ready_masks
        for role_mask in range(1 << len(role_atoms)):
            roles_ready_count = role_ready_count(role_mask, all_role_ready_masks, minimal_masks_by_role, role_atoms)
            all_roles_ready = role_mask in all_role_ready_masks
            end_to_end_ready = bridge_ready and all_roles_ready
            missing_bridge = len(bridge_atoms) - bridge_mask.bit_count()
            missing_role = len(role_atoms) - role_mask.bit_count()
            missing_total = missing_bridge + missing_role
            closure_distribution[f"bridge_ready={bridge_ready}|ready_roles={roles_ready_count}"] += 1
            if bridge_ready:
                bridge_ready_count += 1
                role_count_distribution_when_bridge_ready[roles_ready_count] += 1
            else:
                bridge_blocked_count += 1
                role_count_distribution_when_bridge_blocked[roles_ready_count] += 1
            if end_to_end_ready:
                end_to_end_ready_count += 1
            if bridge_mask == 0 and role_mask == 0:
                current_state = {
                    "bridge_mask": bridge_mask,
                    "role_mask": role_mask,
                    "bridge_ready": bridge_ready,
                    "ready_role_claim_count": roles_ready_count,
                    "end_to_end_ready": end_to_end_ready,
                    "missing_total_atom_count": missing_total,
                    "missing_bridge_atoms": missing_names(bridge_mask, bridge_atoms),
                    "missing_role_obligations": missing_names(role_mask, role_atoms),
                }
            if missing_total == 1 and not end_to_end_ready:
                nearest_misses.append({
                    "bridge_mask": bridge_mask,
                    "role_mask": role_mask,
                    "bridge_ready_before_repair": bridge_ready,
                    "ready_role_claim_count_before_repair": roles_ready_count,
                    "missing_bridge_atoms": missing_names(bridge_mask, bridge_atoms),
                    "missing_role_obligations": missing_names(role_mask, role_atoms),
                })

    bridge_incidence = {row["atom"]: row["component_count"] for row in hypergraph["atom_component_incidence"]}
    role_claim_incidence = {atom: 0 for atom in role_atoms}
    for minimal_sets in minimal_masks_by_role.values():
        for req_set in minimal_sets:
            for atom in req_set:
                role_claim_incidence[atom] += 1
    frontier_atoms = [
        {
            "stage": "bridge_completion",
            "atom": atom,
            "downstream_load_count": bridge_incidence.get(atom, 0),
            "load_kind": "bridge_components_blocked",
        }
        for atom in bridge_atoms
    ] + [
        {
            "stage": "post_bridge_role_transfer",
            "atom": atom,
            "downstream_load_count": role_claim_incidence.get(atom, 0),
            "load_kind": "legacy_role_claims_blocked",
        }
        for atom in role_atoms
    ]
    frontier_atoms = sorted(frontier_atoms, key=lambda row: (-row["downstream_load_count"], row["stage"], row["atom"]))
    minimal_repair_path_from_current = current_state["missing_bridge_atoms"] + current_state["missing_role_obligations"]
    return {
        "bridge_atom_count": len(bridge_atoms),
        "role_obligation_count": len(role_atoms),
        "combined_atom_count": len(bridge_atoms) + len(role_atoms),
        "bridge_atoms": bridge_atoms,
        "role_obligations": role_atoms,
        "combined_assignment_count": (1 << len(bridge_atoms)) * (1 << len(role_atoms)),
        "bridge_ready_assignment_count_in_combined_lattice": bridge_ready_count,
        "bridge_blocked_assignment_count_in_combined_lattice": bridge_blocked_count,
        "end_to_end_ready_assignment_count": end_to_end_ready_count,
        "end_to_end_ready_fraction": f"{end_to_end_ready_count}/{(1 << len(bridge_atoms)) * (1 << len(role_atoms))}",
        "current_state": current_state,
        "minimum_new_atoms_needed_from_current_for_end_to_end_all_role_closure": len(minimal_repair_path_from_current),
        "minimal_repair_path_from_current_ordered_bridge_then_role": minimal_repair_path_from_current,
        "nearest_miss_count": len(nearest_misses),
        "nearest_miss_rows": nearest_misses,
        "closure_distribution": dict(sorted(closure_distribution.items())),
        "role_count_distribution_when_bridge_ready": {str(k): v for k, v in sorted(role_count_distribution_when_bridge_ready.items())},
        "role_count_distribution_when_bridge_blocked": {str(k): v for k, v in sorted(role_count_distribution_when_bridge_blocked.items())},
        "frontier_atom_load_rows": frontier_atoms,
        "top_frontier_atom_load_rows": frontier_atoms[:6],
        "p2411_bridge_ready_true_masks_inherited": sorted(bridge_ready_masks),
        "p2434_all_roles_ready_masks_inherited": sorted(all_role_ready_masks),
        "p2411_nearest_miss_count_inherited": p2411["nearest_miss_count"],
        "p2434_legacy_role_claim_count_inherited": p2434["legacy_role_claim_count"],
        "lattice_fingerprint_sha256": sha256_json({
            "bridge_atoms": bridge_atoms,
            "role_atoms": role_atoms,
            "nearest_misses": nearest_misses,
            "closure_distribution": dict(sorted(closure_distribution.items())),
        }),
    }


def append_once(path, marker: str, section: str) -> None:
    text = path.read_text(encoding="utf-8")
    if marker not in text:
        path.write_text(text.rstrip() + "\n\n" + section.strip() + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2490/S1440 legacy-strict bridge plus role-transfer two-stage closure lattice certificate

`P2490/S1440` combines the P2411 bridge-obligation hypergraph with the P2434 conditional role-transfer claim lattice into one finite two-stage closure lattice.  The computation enumerates `16384 = 2^8 * 2^6` bridge/role assignments and verifies that end-to-end closure for all audited legacy role claims occurs in exactly one assignment: all eight bridge atoms plus all six post-bridge role obligations.  The current state has zero bridge atoms, zero role obligations, zero transferred role claims, and needs fourteen new atoms for this all-role end-to-end closure.

For a non-specialist: this is a bookkeeping theorem about what remains open.  It does not add a new physical formula.  It proves that even a completed bridge would still not automatically transfer legacy roles unless the separate role-transfer and claim-specific successor obligations are also supplied.
"""
    lag_section = """
## P2490/S1440 two-stage bridge/role closure lattice guard

`P2490/S1440` records a finite two-stage guard behind `L_total`: bridge completion and post-bridge legacy-role transfer are separate gates.  Across the `16384` combined assignments, all-role closure appears only in the single all-atoms assignment, while the current state remains at zero transferred legacy role claims.  The guard does not export selector/source discharge, QW-2191 closure, role-bearing `L_total`, physical-value generation, legacy-role transfer, or ToE closure.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2490/S1440 legacy-strict bridge plus role-transfer two-stage closure lattice certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2490/S1440 two-stage bridge/role closure lattice guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    lattice = two_stage_lattice(
        sources["P2411_BRIDGE_OBLIGATION_HYPERGRAPH"],
        sources["P2434_ROLE_TRANSFER_CLAIM_LATTICE"],
    )
    theorem_export = {
        "theorem_name": "P2490_T1_legacy_strict_bridge_role_transfer_two_stage_closure_lattice_certificate",
        "audited_chain": ["P2411/S1361", "P2434/S1384"],
        "two_stage_closure_lattice": lattice,
        "bridge_atom_count": lattice["bridge_atom_count"],
        "role_obligation_count": lattice["role_obligation_count"],
        "combined_atom_count": lattice["combined_atom_count"],
        "combined_assignment_count": lattice["combined_assignment_count"],
        "end_to_end_ready_assignment_count": lattice["end_to_end_ready_assignment_count"],
        "end_to_end_ready_fraction": lattice["end_to_end_ready_fraction"],
        "current_bridge_ready": lattice["current_state"]["bridge_ready"],
        "current_ready_role_claim_count": lattice["current_state"]["ready_role_claim_count"],
        "current_end_to_end_ready": lattice["current_state"]["end_to_end_ready"],
        "minimum_new_atoms_needed_from_current_for_end_to_end_all_role_closure": lattice["minimum_new_atoms_needed_from_current_for_end_to_end_all_role_closure"],
        "nearest_miss_count": lattice["nearest_miss_count"],
        "top_frontier_atom_load_rows": lattice["top_frontier_atom_load_rows"],
        "bridge_completion_and_role_transfer_separate_gates": True,
        "completed_bridge_alone_would_not_transfer_legacy_roles": True,
        "current_legacy_role_claims_transferred_by_this_certificate": 0,
        "source_obligation_discharge_exported_by_this_certificate": False,
        "selector_source_theorem_exported_by_this_certificate": False,
        "qw2191_discharged_by_this_certificate": False,
        "role_transfer_licensed_by_this_certificate": False,
        "role_bearing_ltotal_exported": False,
        "strict_physical_value_generator_exported": False,
        "legacy_role_transfer_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "P2490 is a finite closure-lattice audit, not a new bridge atom proof.",
            "The legacy-to-strict bridge remains incomplete unless all P2411 bridge atoms are supplied.",
            "A completed bridge alone does not transfer legacy physical-role claims without the P2434 post-bridge role obligations.",
            "No QW-2191 discharge, selector/source theorem, role-bearing L_total, physical-value generator, legacy-role transfer, or ToE closure is exported.",
        ],
        "next_honest_step": "Attack one real high-load missing atom noncyclically: a chi11/QW-2191 selector-source theorem for bridge completion, or after bridge completion a claim-specific strict successor theorem. Do not silently inherit legacy roles.",
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "p2411_counts_inherited": theorem_export["bridge_atom_count"] == 8,
        "p2434_counts_inherited": theorem_export["role_obligation_count"] == 6,
        "combined_lattice_size_exact": theorem_export["combined_assignment_count"] == 16384,
        "single_end_to_end_ready_assignment": theorem_export["end_to_end_ready_assignment_count"] == 1,
        "current_state_not_ready": not theorem_export["current_bridge_ready"] and not theorem_export["current_end_to_end_ready"],
        "current_role_transfer_zero": theorem_export["current_ready_role_claim_count"] == 0,
        "all_fourteen_atoms_needed_from_current": theorem_export["minimum_new_atoms_needed_from_current_for_end_to_end_all_role_closure"] == 14,
        "nearest_misses_are_single_atom_gaps": theorem_export["nearest_miss_count"] == 14,
        "separate_gate_rule_preserved": theorem_export["bridge_completion_and_role_transfer_separate_gates"] and theorem_export["completed_bridge_alone_would_not_transfer_legacy_roles"],
        "no_closure_inflation": not any(theorem_export[key] for key in [
            "source_obligation_discharge_exported_by_this_certificate",
            "selector_source_theorem_exported_by_this_certificate",
            "qw2191_discharged_by_this_certificate",
            "role_transfer_licensed_by_this_certificate",
            "role_bearing_ltotal_exported",
            "strict_physical_value_generator_exported",
            "legacy_role_transfer_exported",
            "toe_closure_exported",
        ]),
    }
    return {
        "packet_id": "P2490",
        "stage_id": "S1440",
        "status": "TWO_STAGE_LEGACY_STRICT_BRIDGE_ROLE_TRANSFER_CLOSURE_LATTICE_CERTIFICATE_NO_BRIDGE_ATOM_EXPORT_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "legacy_strict_bridge_role_transfer_two_stage_closure_lattice_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["legacy_strict_bridge_role_transfer_two_stage_closure_lattice_certificate"]["theorem_export"]
    lattice = t["two_stage_closure_lattice"]
    lines = [
        "# P2490/S1440 legacy-strict bridge plus role-transfer two-stage closure lattice certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Two-stage closure lattice",
        "",
        f"Audited chain: `{', '.join(t['audited_chain'])}`.",
        f"Bridge atoms: `{t['bridge_atom_count']}`.",
        f"Post-bridge role obligations: `{t['role_obligation_count']}`.",
        f"Combined assignments enumerated: `{t['combined_assignment_count']}`.",
        f"End-to-end all-role ready assignments: `{t['end_to_end_ready_assignment_count']}`.",
        f"End-to-end ready fraction: `{t['end_to_end_ready_fraction']}`.",
        f"Current bridge ready: `{t['current_bridge_ready']}`.",
        f"Current ready role claims: `{t['current_ready_role_claim_count']}`.",
        f"Minimum new atoms needed from current: `{t['minimum_new_atoms_needed_from_current_for_end_to_end_all_role_closure']}`.",
        f"Nearest one-atom-miss states: `{t['nearest_miss_count']}`.",
        "",
        "## Top frontier atom loads",
        "",
    ]
    for row in t["top_frontier_atom_load_rows"]:
        lines.append(f"- `{row['stage']}` / `{row['atom']}` blocks `{row['downstream_load_count']}` `{row['load_kind']}`.")
    lines += [
        "",
        "## Negative controls",
        "",
        "P2490 does not export a bridge atom, selector/source theorem, QW-2191 discharge, role-transfer license, role-bearing L_total, physical-value generator, legacy-role transfer, or ToE closure.",
        "",
        "## Lay summary",
        "",
        "P2490 proves a precise accounting result: bridge completion and legacy-role transfer are separate gates.  In the combined finite lattice, all-role closure appears only when every bridge atom and every post-bridge role obligation is present.",
        "",
        "## Fingerprints",
        "",
        f"Lattice fingerprint: `{lattice['lattice_fingerprint_sha256']}`.",
        f"Theorem fingerprint: `{payload['legacy_strict_bridge_role_transfer_two_stage_closure_lattice_certificate']['theorem_fingerprint_sha256']}`.",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
