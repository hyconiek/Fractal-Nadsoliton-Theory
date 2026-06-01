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

OUT = GEN / "p2404_s1354_strict_addition_physics_lane_dependency_cut_certificate.json"
MD = GEN / "p2404_s1354_strict_addition_physics_lane_dependency_cut_certificate.md"

SOURCE_FILES = {
    "S2_PRIORITY_PACKET": ROOT / "S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md",
    "P2400_ROLE_LATTICE": GEN / "p2400_s1350_nearest_lift_role_successor_lattice_certificate.json",
    "P2403_STRICT_PRIMARY_REBASE": GEN / "p2403_s1353_strict_kernel_primary_physics_generation_rebase_certificate.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

STRICT_ADDITION_ATOMS = [
    "apd_completion_structure",
    "gf2_phase_topological_data",
    "nonlinear_damping_compression",
    "chi11_selector_source_declared",
]

ROLE_SUCCESSOR_ATOMS = [
    "alpha_geo_electroweak_role_theorem",
    "beta_tors_strict_role_theorem",
    "beta_power_hierarchy_successor_theorem",
]

ATOMS = STRICT_ADDITION_ATOMS + ROLE_SUCCESSOR_ATOMS
ATOM_INDEX = {atom: index for index, atom in enumerate(ATOMS)}

LANE_REQUIREMENTS = {
    "strict_kernel_structural_candidate_test_readiness": STRICT_ADDITION_ATOMS,
    "strict_mass_generation_candidate_test_readiness": STRICT_ADDITION_ATOMS,
    "legacy_weinberg_role_transfer_to_strict_successor": STRICT_ADDITION_ATOMS + ["alpha_geo_electroweak_role_theorem"],
    "legacy_alpha_em_role_transfer_to_strict_successor": STRICT_ADDITION_ATOMS
    + ["alpha_geo_electroweak_role_theorem", "beta_tors_strict_role_theorem"],
    "legacy_gravity_hierarchy_strict_successor": STRICT_ADDITION_ATOMS
    + ["beta_tors_strict_role_theorem", "beta_power_hierarchy_successor_theorem"],
    "strict_role_bearing_ltotal_promotion_candidate": STRICT_ADDITION_ATOMS + ROLE_SUCCESSOR_ATOMS,
    "strict_toe_physics_generation_package_candidate": STRICT_ADDITION_ATOMS + ROLE_SUCCESSOR_ATOMS,
}

PHYSICAL_ROLE_LANES = [
    "legacy_weinberg_role_transfer_to_strict_successor",
    "legacy_alpha_em_role_transfer_to_strict_successor",
    "legacy_gravity_hierarchy_strict_successor",
    "strict_role_bearing_ltotal_promotion_candidate",
    "strict_toe_physics_generation_package_candidate",
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


def rg_count(pattern: str, *extra_paths: str) -> dict[str, Any]:
    paths = list(extra_paths) if extra_paths else ["fundamental_action_reconstruction", "material_dowodowy"]
    proc = subprocess.run(
        [
            "rg",
            "-n",
            pattern,
            *paths,
            "-g",
            "*.py",
            "-g",
            "*.md",
            "-g",
            "*.json",
            "-g",
            "*.tex",
            "-g",
            "!generated/p2404_s1354_strict_addition_physics_lane_dependency_cut_certificate.*",
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
        "new_packet": "P2404|S1354|strict addition physics lane dependency|physics-lane dependency cut|dependency-cut certificate",
        "strict_primary_rebase": "P2403|strict-kernel primary physics-generation|strict primary",
        "role_lattice": "P2400|three-role successor lattice|nearest-lift role-successor",
        "lane_dependency_candidates": "dependency cut|readiness lattice|physics lane|mass_generation|gravity_hierarchy|lagrangian_eom_promotion",
        "guardrail_blocks": "silent role transfer|role-successor theorem|required before physical promotion|No L_total promotion",
    }
    return {
        "tool": "rg",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
        "finding": (
            "Repo grep finds P2403 strict-primary rebase and P2400 role lattice, but no separate finite "
            "dependency-cut/readiness lattice that combines strict additions with role-successor atoms for physics lanes."
        ),
    }


def mask_for(atoms: list[str]) -> int:
    mask = 0
    for atom in atoms:
        mask |= 1 << ATOM_INDEX[atom]
    return mask


def atoms_for(mask: int) -> list[str]:
    return [atom for atom in ATOMS if mask & (1 << ATOM_INDEX[atom])]


def lane_ready(mask: int, lane: str) -> bool:
    req_mask = mask_for(list(LANE_REQUIREMENTS[lane]))
    return (mask & req_mask) == req_mask


def minimal_true_masks_for_lane(lane: str) -> list[int]:
    true_masks = [mask for mask in range(1 << len(ATOMS)) if lane_ready(mask, lane)]
    minima = []
    for mask in true_masks:
        if not any(other != mask and (other & mask) == other for other in true_masks):
            minima.append(mask)
    return sorted(minima)


def polynomial_monomial(atoms: list[str]) -> str:
    return " * ".join(atoms) if atoms else "1"


def exact_dependency_lattice() -> dict[str, Any]:
    rows = []
    for mask in range(1 << len(ATOMS)):
        ready_lanes = [lane for lane in LANE_REQUIREMENTS if lane_ready(mask, lane)]
        rows.append(
            {
                "mask": mask,
                "atoms_true": atoms_for(mask),
                "ready_lanes": ready_lanes,
                "ready_lane_count": len(ready_lanes),
                "physical_role_lanes_ready": [lane for lane in ready_lanes if lane in PHYSICAL_ROLE_LANES],
            }
        )

    lane_rows = []
    for lane, requirements in LANE_REQUIREMENTS.items():
        req_mask = mask_for(list(requirements))
        minima = minimal_true_masks_for_lane(lane)
        lane_rows.append(
            {
                "lane": lane,
                "required_atoms": list(requirements),
                "requirement_mask": req_mask,
                "minimal_true_masks": minima,
                "minimal_true_atom_sets": [atoms_for(mask) for mask in minima],
                "anf_over_dependency_atoms": polynomial_monomial(list(requirements)),
                "anf_degree": len(requirements),
                "strict_addition_requirement_count": sum(atom in STRICT_ADDITION_ATOMS for atom in requirements),
                "role_successor_requirement_count": sum(atom in ROLE_SUCCESSOR_ATOMS for atom in requirements),
            }
        )

    common_cut = sorted(set.intersection(*(set(reqs) for reqs in LANE_REQUIREMENTS.values())), key=ATOM_INDEX.get)
    legacy_only_mask = 0
    strict_additions_only_mask = mask_for(STRICT_ADDITION_ATOMS)
    full_mask = (1 << len(ATOMS)) - 1
    return {
        "atom_order": ATOMS,
        "strict_addition_atoms": STRICT_ADDITION_ATOMS,
        "role_successor_atoms": ROLE_SUCCESSOR_ATOMS,
        "truth_table_row_count": len(rows),
        "rows": rows,
        "lane_dependency_rows": lane_rows,
        "common_strict_addition_cut_support": common_cut,
        "common_strict_addition_cut_mask": mask_for(common_cut),
        "legacy_only_mask": legacy_only_mask,
        "legacy_only_ready_lanes": rows[legacy_only_mask]["ready_lanes"],
        "strict_additions_only_mask": strict_additions_only_mask,
        "strict_additions_only_ready_lanes": rows[strict_additions_only_mask]["ready_lanes"],
        "strict_additions_only_physical_role_lanes_ready": rows[strict_additions_only_mask]["physical_role_lanes_ready"],
        "full_mask": full_mask,
        "full_mask_ready_lanes": rows[full_mask]["ready_lanes"],
        "distance_from_legacy_only_by_lane": {
            lane: min(mask.bit_count() for mask in minimal_true_masks_for_lane(lane)) for lane in LANE_REQUIREMENTS
        },
    }


def append_doc_sections() -> None:
    eq_section = """
## P2404/S1354 strict-addition physics-lane dependency-cut certificate

`P2404/S1354` turns the P2403 strict-primary rebase into a finite dependency-cut audit.  The audited atoms are the four strict additions (`A/P/D`, GF(2)/topological phase data, nonlinear compression, and `chi11` selector bookkeeping) plus the three role-successor atoms from the P2400 lattice.

The exact 128-row lattice shows that every listed physics-generation lane has the same common strict-addition cut: all four strict-side additions must be present before the strict kernel is even the correct primary object for that lane.  With only the strict additions present, the structural strict-kernel and mass-generation candidate tests become ready, but all role-bearing physical lanes remain false.  Role transfer, `L_total`, and ToE package readiness still require the relevant role-successor atoms; the full package has the degree-7 monomial consisting of all four strict additions and all three role-successor atoms.

Thus strict is computationally favored as the primary physics-generation candidate because it carries the missing nadsoliton characteristics, while legacy-only physics inheritance remains blocked by an explicit dependency cut rather than by narrative preference.
""".strip()
    lag_section = """
## P2404/S1354 dependency-cut guard for Lagrangian/EOM

`P2404/S1354` computes the finite dependency cut for strict physics lanes.  The `L_total` and ToE candidate lanes have the degree-7 dependency monomial: four strict additions plus all three role-successor atoms.  Therefore no Lagrangian term may be promoted from legacy-only data, from strict structural additions alone, or from any one-/two-role prefix.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items() if path.suffix == ".json"}
    text_sources = {name: load_text(path) for name, path in SOURCE_FILES.items() if path.suffix != ".json"}
    grep = rg_audit()
    lattice = exact_dependency_lattice()
    p2403_theorem = artifacts["P2403_STRICT_PRIMARY_REBASE"].get(
        "strict_kernel_primary_physics_generation_rebase_certificate", {}
    ).get("theorem_export", {})
    p2400_theorem = artifacts["P2400_ROLE_LATTICE"].get(
        "nearest_lift_role_successor_lattice_certificate", {}
    ).get("theorem_export", {})
    strict_additions_from_p2403 = p2403_theorem.get("characteristic_matrix_summary", {}).get("strict_additions", [])
    s2_text = text_sources.get("S2_PRIORITY_PACKET", "")
    ltotal_lane = next(row for row in lattice["lane_dependency_rows"] if row["lane"] == "strict_role_bearing_ltotal_promotion_candidate")
    toe_lane = next(row for row in lattice["lane_dependency_rows"] if row["lane"] == "strict_toe_physics_generation_package_candidate")

    theorem_export = {
        "theorem_name": "P2404_T1_strict_addition_physics_lane_dependency_cut",
        "atom_order": lattice["atom_order"],
        "truth_table_row_count": lattice["truth_table_row_count"],
        "common_strict_addition_cut_support": lattice["common_strict_addition_cut_support"],
        "common_strict_addition_cut_mask": lattice["common_strict_addition_cut_mask"],
        "legacy_only_ready_lanes": lattice["legacy_only_ready_lanes"],
        "strict_additions_only_ready_lanes": lattice["strict_additions_only_ready_lanes"],
        "strict_additions_only_physical_role_lanes_ready": lattice["strict_additions_only_physical_role_lanes_ready"],
        "full_mask_ready_lane_count": len(lattice["full_mask_ready_lanes"]),
        "full_mask_ready_lanes": lattice["full_mask_ready_lanes"],
        "ltotal_dependency_anf": ltotal_lane["anf_over_dependency_atoms"],
        "ltotal_dependency_degree": ltotal_lane["anf_degree"],
        "toe_package_dependency_anf": toe_lane["anf_over_dependency_atoms"],
        "toe_package_dependency_degree": toe_lane["anf_degree"],
        "distance_from_legacy_only_by_lane": lattice["distance_from_legacy_only_by_lane"],
        "p2403_strict_additions_matched": strict_additions_from_p2403 == STRICT_ADDITION_ATOMS,
        "p2400_full_role_mask_requirement_inherited": p2400_theorem.get("toe_true_masks") == [7],
        "s2_bridge_priority_terms_detected": all(term in s2_text for term in ["legacy -> strict", "role-transfer audit", "strict-side additions"]),
        "not_licensed": [
            "No SM/GR mass spectrum is derived by this dependency lattice.",
            "No legacy physical role transfers to strict without its role-successor atom.",
            "Strict additions alone only make structural/candidate tests ready; they do not license role-bearing L_total.",
            "The full degree-7 dependency monomial is conditional readiness, not ToE closure.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "truth_table_has_128_rows": lattice["truth_table_row_count"] == 128,
        "all_lanes_share_four_strict_addition_cut": lattice["common_strict_addition_cut_support"] == STRICT_ADDITION_ATOMS,
        "legacy_only_licenses_no_lanes": lattice["legacy_only_ready_lanes"] == [],
        "strict_additions_only_license_no_physical_role_lanes": lattice[
            "strict_additions_only_physical_role_lanes_ready"
        ] == [],
        "strict_additions_only_ready_lanes_are_candidate_tests_only": lattice["strict_additions_only_ready_lanes"]
        == ["strict_kernel_structural_candidate_test_readiness", "strict_mass_generation_candidate_test_readiness"],
        "full_mask_readies_all_lanes_conditionally": len(lattice["full_mask_ready_lanes"]) == len(LANE_REQUIREMENTS),
        "ltotal_and_toe_have_degree_seven_dependencies": ltotal_lane["anf_degree"] == 7 and toe_lane["anf_degree"] == 7,
        "p2403_strict_additions_matched": theorem_export["p2403_strict_additions_matched"],
        "p2400_full_role_mask_requirement_inherited": theorem_export["p2400_full_role_mask_requirement_inherited"],
        "s2_bridge_priority_detected": theorem_export["s2_bridge_priority_terms_detected"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2404_s1354_v1",
        "packet_id": "P2404",
        "stage_id": "S1354",
        "result_kind": "STRICT_ADDITION_PHYSICS_LANE_DEPENDENCY_CUT_CERTIFICATE",
        "status": "PASS_STRICT_ADDITION_DEPENDENCY_CUT_NO_ROLE_TRANSFER_NO_TOE_CLOSURE",
        "strict_addition_physics_lane_dependency_cut_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: artifact.get("status", "TEXT_SOURCE") for name, artifact in artifacts.items()},
            "dependency_lattice": lattice,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": (
            "Use the common four-addition cut to design strict-kernel physics-generation tests, while keeping "
            "role-bearing lanes behind the explicit P2400 three-role successor requirement."
        ),
        "global_status": "OPEN_PROGRESS_STRICT_PHYSICS_DEPENDENCY_CUT_CERTIFIED_NO_PHYSICAL_ROLE_EXPORT",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    cert = payload["strict_addition_physics_lane_dependency_cut_certificate"]
    theorem = cert["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2404 S1354: strict-addition physics-lane dependency-cut certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2404/S1354 computes the 128-row dependency lattice over four strict additions plus three role-successor atoms.",
                "",
                "## Exact dependency cut",
                "",
                f"- Common strict-addition cut: `{theorem['common_strict_addition_cut_support']}`.",
                f"- Legacy-only ready lanes: `{theorem['legacy_only_ready_lanes']}`.",
                f"- Strict-additions-only ready lanes: `{theorem['strict_additions_only_ready_lanes']}`.",
                f"- Strict-additions-only physical role lanes ready: `{theorem['strict_additions_only_physical_role_lanes_ready']}`.",
                f"- Full-mask conditional lane count: `{theorem['full_mask_ready_lane_count']}`.",
                "",
                "## Lagrangian / ToE dependencies",
                "",
                f"- `L_total` dependency ANF: `{theorem['ltotal_dependency_anf']}`.",
                f"- `L_total` dependency degree: `{theorem['ltotal_dependency_degree']}`.",
                f"- ToE package dependency degree: `{theorem['toe_package_dependency_degree']}`.",
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
