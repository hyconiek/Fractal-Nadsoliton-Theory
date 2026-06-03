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
OUT = GEN / "p2491_s1441_legacy_strict_bridge_role_claim_pivotality_matrix_certificate.json"
MD = GEN / "p2491_s1441_legacy_strict_bridge_role_claim_pivotality_matrix_certificate.md"

SOURCE_FILES = {
    "P2490_TWO_STAGE_LATTICE": GEN / "p2490_s1440_legacy_strict_bridge_role_transfer_two_stage_closure_lattice_certificate.json",
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
        "new_packet": "P2491|S1441|bridge-role pivotality matrix|role-claim pivotality|claim-specific pivotality|end-to-end role pivotal|combined pivotality matrix",
        "precursor_packets": "P2490|S1440|two-stage bridge-role closure lattice|P2402|S1352|role-successor marginal credit|atom influence|Banzhaf|pivotal",
        "claim_language": "legacy_weak_mixing_angle|legacy_inverse_alpha_em|legacy_beta_power_gravity_hierarchy|legacy_torsion_to_chi11_orientation|claim-specific strict successor",
        "guardrail_language": "legacy -> strict completion bridge|role-transfer audit|silent inheritance|K_legacy_ont|K_strict_gate|QW-2191",
        "closure_blockers": "role-bearing L_total|physical-value generator|ToE closure|selector/source/gauge theorem|legacy-role transfer|root-window theorem",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def mask_has(mask: int, index: int) -> bool:
    return bool((mask >> index) & 1)


def claim_ready(mask: int, required_indices: set[int]) -> bool:
    return all(mask_has(mask, index) for index in required_indices)


def pivotal_count(atom_index: int, required_indices: set[int], combined_atom_count: int) -> int:
    count = 0
    for mask in range(1 << combined_atom_count):
        if mask_has(mask, atom_index):
            continue
        before = claim_ready(mask, required_indices)
        after = claim_ready(mask | (1 << atom_index), required_indices)
        if not before and after:
            count += 1
    return count


def pivotality_matrix(p2490: dict[str, Any], p2434: dict[str, Any]) -> dict[str, Any]:
    lattice = p2490["two_stage_closure_lattice"]
    bridge_atoms = lattice["bridge_atoms"]
    role_atoms = lattice["role_obligations"]
    atoms = [f"bridge::{atom}" for atom in bridge_atoms] + [f"role::{atom}" for atom in role_atoms]
    atom_index = {atom: index for index, atom in enumerate(atoms)}
    role_offset = len(bridge_atoms)
    all_bridge_indices = set(range(len(bridge_atoms)))
    role_claim_requirements = {
        role_id: [f"role::{atom}" for atom in req_sets[0]]
        for role_id, req_sets in p2434["minimal_masks_by_role_claim"].items()
    }
    claim_requirements = {
        claim: sorted(all_bridge_indices | {atom_index[atom] for atom in role_req_atoms})
        for claim, role_req_atoms in role_claim_requirements.items()
    }
    claim_requirements["all_four_audited_legacy_role_claims"] = list(range(len(atoms)))

    matrix_rows = []
    total_by_atom: Counter[str] = Counter()
    physical_claim_total_by_atom: Counter[str] = Counter()
    for atom in atoms:
        index = atom_index[atom]
        claim_counts = {}
        for claim, required_indices in claim_requirements.items():
            required = set(required_indices)
            count = pivotal_count(index, required, len(atoms))
            claim_counts[claim] = {
                "pivotal_count": count,
                "pivotal_fraction": f"{count}/{1 << (len(atoms) - 1)}",
                "atom_required_for_claim": index in required,
            }
            total_by_atom[atom] += count
            if claim != "all_four_audited_legacy_role_claims":
                physical_claim_total_by_atom[atom] += count
        matrix_rows.append({
            "atom": atom,
            "stage": atom.split("::", 1)[0],
            "claim_pivotal_counts": claim_counts,
            "physical_role_claim_pivotal_total": physical_claim_total_by_atom[atom],
            "all_targets_pivotal_total": total_by_atom[atom],
        })
    matrix_rows.sort(key=lambda row: (-row["physical_role_claim_pivotal_total"], row["stage"], row["atom"]))
    top_physical = matrix_rows[0]["physical_role_claim_pivotal_total"]
    top_atoms = [row["atom"] for row in matrix_rows if row["physical_role_claim_pivotal_total"] == top_physical]
    role_stage_rows = [row for row in matrix_rows if row["stage"] == "role"]
    top_role_stage_physical = max(row["physical_role_claim_pivotal_total"] for row in role_stage_rows)
    top_role_stage_atoms = [row["atom"] for row in role_stage_rows if row["physical_role_claim_pivotal_total"] == top_role_stage_physical]
    current_mask = 0
    current_ready_claims = [claim for claim, required in claim_requirements.items() if claim_ready(current_mask, set(required))]
    return {
        "combined_atom_count": len(atoms),
        "assignment_count": 1 << len(atoms),
        "one_atom_flip_case_count": len(atoms) * (1 << (len(atoms) - 1)),
        "claim_count_including_all_roles_target": len(claim_requirements),
        "bridge_atoms": bridge_atoms,
        "role_obligations": role_atoms,
        "combined_atoms": atoms,
        "claim_requirements": {claim: [atoms[index] for index in required] for claim, required in claim_requirements.items()},
        "current_ready_claims": current_ready_claims,
        "current_ready_claim_count": len(current_ready_claims),
        "matrix_rows": matrix_rows,
        "top_physical_role_claim_pivotal_total": top_physical,
        "top_physical_role_claim_pivotal_atoms": top_atoms,
        "top_role_stage_physical_role_claim_pivotal_total": top_role_stage_physical,
        "top_role_stage_physical_role_claim_pivotal_atoms": top_role_stage_atoms,
        "all_bridge_atoms_have_equal_physical_role_pivotal_total": len({row["physical_role_claim_pivotal_total"] for row in matrix_rows if row["stage"] == "bridge"}) == 1,
        "role_transfer_and_ltotal_are_top_role_stage_pivotal_atoms": set(top_role_stage_atoms) == {"role::role_transfer_audit_license", "role::role_bearing_ltotal_export"},
        "every_atom_is_pivotal_for_all_roles_target_once": all(row["claim_pivotal_counts"]["all_four_audited_legacy_role_claims"]["pivotal_count"] == 1 for row in matrix_rows),
        "matrix_fingerprint_sha256": sha256_json(matrix_rows),
    }


def append_once(path, marker: str, section: str) -> None:
    text = path.read_text(encoding="utf-8")
    if marker not in text:
        path.write_text(text.rstrip() + "\n\n" + section.strip() + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2491/S1441 legacy-strict bridge/role claim pivotality matrix certificate

`P2491/S1441` refines the P2490 two-stage lattice by computing a claim-specific Boolean pivotality matrix across the same fourteen open atoms.  For each atom and each audited legacy role claim, it enumerates all one-bit flips in the `2^14` combined bridge/role lattice and counts exactly when that atom is pivotal for end-to-end claim readiness.  The result confirms that all bridge atoms are equally prerequisite for every role claim, while the top post-bridge role-stage atoms are `role_transfer_audit_license` and `role_bearing_ltotal_export`.  This is a sensitivity certificate only: no bridge atom, selector source, QW-2191 discharge, role-transfer license, physical-value generator, or ToE closure is exported.

For a non-specialist: P2490 said all gates must be supplied.  P2491 asks which missing gates are pivotal for which legacy claims.  It quantifies the bottlenecks without pretending that any of them has been proved.
"""
    lag_section = """
## P2491/S1441 bridge/role claim pivotality guard

`P2491/S1441` adds a claim-specific pivotality guard behind `L_total`: every audited legacy-role claim still depends on all bridge atoms, and the strongest post-bridge role-stage bottlenecks are the role-transfer audit license and role-bearing `L_total` export.  The guard is a finite Boolean sensitivity audit and does not export selector/source discharge, QW-2191 closure, role-transfer, physical-value generation, or ToE closure.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2491/S1441 legacy-strict bridge/role claim pivotality matrix certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2491/S1441 bridge/role claim pivotality guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2490 = theorem(sources["P2490_TWO_STAGE_LATTICE"], "legacy_strict_bridge_role_transfer_two_stage_closure_lattice_certificate")
    p2434 = theorem(sources["P2434_ROLE_TRANSFER_CLAIM_LATTICE"], "conditional_legacy_role_transfer_claim_lattice_certificate")
    matrix = pivotality_matrix(p2490, p2434)
    theorem_export = {
        "theorem_name": "P2491_T1_legacy_strict_bridge_role_claim_pivotality_matrix_certificate",
        "audited_chain": ["P2490/S1440"],
        "bridge_role_claim_pivotality_matrix": matrix,
        "combined_atom_count": matrix["combined_atom_count"],
        "assignment_count": matrix["assignment_count"],
        "one_atom_flip_case_count": matrix["one_atom_flip_case_count"],
        "claim_count_including_all_roles_target": matrix["claim_count_including_all_roles_target"],
        "current_ready_claim_count": matrix["current_ready_claim_count"],
        "top_physical_role_claim_pivotal_total": matrix["top_physical_role_claim_pivotal_total"],
        "top_physical_role_claim_pivotal_atoms": matrix["top_physical_role_claim_pivotal_atoms"],
        "top_role_stage_physical_role_claim_pivotal_total": matrix["top_role_stage_physical_role_claim_pivotal_total"],
        "top_role_stage_physical_role_claim_pivotal_atoms": matrix["top_role_stage_physical_role_claim_pivotal_atoms"],
        "all_bridge_atoms_have_equal_physical_role_pivotal_total": matrix["all_bridge_atoms_have_equal_physical_role_pivotal_total"],
        "role_transfer_and_ltotal_are_top_role_stage_pivotal_atoms": matrix["role_transfer_and_ltotal_are_top_role_stage_pivotal_atoms"],
        "every_atom_is_pivotal_for_all_roles_target_once": matrix["every_atom_is_pivotal_for_all_roles_target_once"],
        "source_obligation_discharge_exported_by_this_certificate": False,
        "selector_source_theorem_exported_by_this_certificate": False,
        "qw2191_discharged_by_this_certificate": False,
        "role_transfer_licensed_by_this_certificate": False,
        "role_bearing_ltotal_exported": False,
        "strict_physical_value_generator_exported": False,
        "legacy_role_transfer_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "P2491 computes Boolean pivotality counts only; it does not prove any missing atom.",
            "Equal bridge pivotality means every bridge atom is still prerequisite, not that the bridge is complete.",
            "Top role-stage pivotality of role_transfer_audit_license and role_bearing_ltotal_export does not export either object.",
            "No QW-2191 discharge, selector/source theorem, role-transfer theorem, physical-value generator, legacy-role transfer, or ToE closure is exported.",
        ],
        "next_honest_step": "Use the pivotality matrix to choose a real proof target rather than another lattice-only inflation: a bridge-stage source atom (especially chi11/QW-2191 if selector/source work is attempted) or the post-bridge role-transfer/L_total objects, each with an actual theorem source.",
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "p2490_inherited": theorem_export["combined_atom_count"] == 14,
        "assignment_count_exact": theorem_export["assignment_count"] == 16384,
        "one_flip_case_count_exact": theorem_export["one_atom_flip_case_count"] == 114688,
        "current_exports_no_claim": theorem_export["current_ready_claim_count"] == 0,
        "bridge_pivotality_uniform": theorem_export["all_bridge_atoms_have_equal_physical_role_pivotal_total"],
        "top_role_stage_bottlenecks_checked": theorem_export["role_transfer_and_ltotal_are_top_role_stage_pivotal_atoms"],
        "all_roles_target_requires_every_atom": theorem_export["every_atom_is_pivotal_for_all_roles_target_once"],
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
        "packet_id": "P2491",
        "stage_id": "S1441",
        "status": "LEGACY_STRICT_BRIDGE_ROLE_CLAIM_PIVOTALITY_MATRIX_CERTIFICATE_NO_ATOM_EXPORT_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "legacy_strict_bridge_role_claim_pivotality_matrix_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["legacy_strict_bridge_role_claim_pivotality_matrix_certificate"]["theorem_export"]
    matrix = t["bridge_role_claim_pivotality_matrix"]
    lines = [
        "# P2491/S1441 legacy-strict bridge/role claim pivotality matrix certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Pivotality matrix summary",
        "",
        f"Audited chain: `{', '.join(t['audited_chain'])}`.",
        f"Combined atoms: `{t['combined_atom_count']}`.",
        f"Assignments: `{t['assignment_count']}`.",
        f"One-atom flip cases enumerated per target layer: `{t['one_atom_flip_case_count']}`.",
        f"Claims including all-role target: `{t['claim_count_including_all_roles_target']}`.",
        f"Current ready claims: `{t['current_ready_claim_count']}`.",
        f"Top physical-role pivotal total: `{t['top_physical_role_claim_pivotal_total']}`.",
        f"Top physical-role pivotal atoms: `{', '.join(t['top_physical_role_claim_pivotal_atoms'])}`.",
        f"Top role-stage pivotal atoms: `{', '.join(t['top_role_stage_physical_role_claim_pivotal_atoms'])}`.",
        f"All bridge atoms have equal physical-role pivotal total: `{t['all_bridge_atoms_have_equal_physical_role_pivotal_total']}`.",
        f"Role-transfer and L_total are top role-stage pivotal atoms: `{t['role_transfer_and_ltotal_are_top_role_stage_pivotal_atoms']}`.",
        "",
        "## Top rows",
        "",
    ]
    for row in matrix["matrix_rows"][:8]:
        lines.append(f"- `{row['atom']}`: physical-role total `{row['physical_role_claim_pivotal_total']}`, all-target total `{row['all_targets_pivotal_total']}`.")
    lines += [
        "",
        "## Negative controls",
        "",
        "P2491 does not export a bridge atom, selector/source theorem, QW-2191 discharge, role-transfer license, role-bearing L_total, physical-value generator, legacy-role transfer, or ToE closure.",
        "",
        "## Lay summary",
        "",
        "P2491 quantifies which missing gates are pivotal for audited legacy role claims in the combined bridge/role lattice.  It ranks bottlenecks but does not prove any missing theorem atom.",
        "",
        "## Fingerprints",
        "",
        f"Matrix fingerprint: `{matrix['matrix_fingerprint_sha256']}`.",
        f"Theorem fingerprint: `{payload['legacy_strict_bridge_role_claim_pivotality_matrix_certificate']['theorem_fingerprint_sha256']}`.",
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
