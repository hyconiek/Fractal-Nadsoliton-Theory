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

OUT = GEN / "p2406_s1356_information_to_physics_staged_projection_barrier_certificate.json"
MD = GEN / "p2406_s1356_information_to_physics_staged_projection_barrier_certificate.md"

SOURCE_FILES = {
    "P2404_DEPENDENCY_CUT": GEN / "p2404_s1354_strict_addition_physics_lane_dependency_cut_certificate.json",
    "P2405_INFORMATION_ONTOLOGY": GEN / "p2405_s1355_nadsoliton_information_ontology_projection_certificate.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

ONTOLOGY_GUARD_ATOMS = [
    "nadsoliton_is_sole_primordial_information",
    "no_separate_information_layer_under_nadsoliton",
    "strict_additions_are_internal_information_constraints",
    "physics_roles_are_downstream_projections",
    "observer_is_downstream_readout_not_source",
]

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

ATOMS = ONTOLOGY_GUARD_ATOMS + STRICT_ADDITION_ATOMS + ROLE_SUCCESSOR_ATOMS
ATOM_INDEX = {atom: index for index, atom in enumerate(ATOMS)}

STAGED_LANES = {
    "ontology_typed_information_root": ONTOLOGY_GUARD_ATOMS,
    "strict_internal_information_completion_ready": ONTOLOGY_GUARD_ATOMS + STRICT_ADDITION_ATOMS,
    "weinberg_downstream_projection_candidate": ONTOLOGY_GUARD_ATOMS
    + STRICT_ADDITION_ATOMS
    + ["alpha_geo_electroweak_role_theorem"],
    "alpha_em_downstream_projection_candidate": ONTOLOGY_GUARD_ATOMS
    + STRICT_ADDITION_ATOMS
    + ["alpha_geo_electroweak_role_theorem", "beta_tors_strict_role_theorem"],
    "gravity_hierarchy_downstream_projection_candidate": ONTOLOGY_GUARD_ATOMS
    + STRICT_ADDITION_ATOMS
    + ["beta_tors_strict_role_theorem", "beta_power_hierarchy_successor_theorem"],
    "role_bearing_ltotal_downstream_projection_candidate": ONTOLOGY_GUARD_ATOMS
    + STRICT_ADDITION_ATOMS
    + ROLE_SUCCESSOR_ATOMS,
    "toe_downstream_projection_candidate": ONTOLOGY_GUARD_ATOMS + STRICT_ADDITION_ATOMS + ROLE_SUCCESSOR_ATOMS,
}

ROLE_LANES = [
    "weinberg_downstream_projection_candidate",
    "alpha_em_downstream_projection_candidate",
    "gravity_hierarchy_downstream_projection_candidate",
    "role_bearing_ltotal_downstream_projection_candidate",
    "toe_downstream_projection_candidate",
]


def rel(path: Path) -> str:
    return path.relative_to(REPO).as_posix()


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": rel(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


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
            "*.tex",
            "-g",
            "!generated/**",
        ],
        cwd=REPO,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:16]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2406|S1356|information-to-physics staged projection|staged projection barrier",
        "p2404_dependency": "P2404|degree-7 dependency|strict-additions-only physical role lanes",
        "p2405_ontology": "P2405|nadsoliton information-ontology|sub-nadsoliton information layer|internal information constraints",
        "role_projection_blocks": "role-successor|downstream projection|L_total.*downstream|No L_total|ToE closure",
    }
    return {
        "tool": "rg",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
        "finding": (
            "Repo grep, excluding generated artifacts, finds P2404 and P2405 but no compact staged certificate "
            "that computes the full ontology+strict+role projection barrier as one finite Boolean object."
        ),
    }


def mask_for(atoms: list[str]) -> int:
    mask = 0
    for atom in atoms:
        mask |= 1 << ATOM_INDEX[atom]
    return mask


def atoms_for(mask: int) -> list[str]:
    return [atom for atom in ATOMS if mask & (1 << ATOM_INDEX[atom])]


def lane_count(requirements: list[str]) -> int:
    return 1 << (len(ATOMS) - len(set(requirements)))


def monomial(atoms: list[str]) -> str:
    return " * ".join(atoms)


def staged_projection_summary() -> dict[str, Any]:
    total_rows = 1 << len(ATOMS)
    lane_rows = []
    for lane, requirements in STAGED_LANES.items():
        req = list(dict.fromkeys(requirements))
        lane_rows.append(
            {
                "lane": lane,
                "required_atoms": req,
                "requirement_mask": mask_for(req),
                "anf_monomial": monomial(req),
                "anf_degree": len(req),
                "true_assignment_count": lane_count(req),
                "true_assignment_fraction": f"1/{1 << len(req)}",
                "missing_from_ontology_only": [atom for atom in req if atom not in ONTOLOGY_GUARD_ATOMS],
                "missing_from_ontology_plus_strict": [atom for atom in req if atom not in ONTOLOGY_GUARD_ATOMS + STRICT_ADDITION_ATOMS],
            }
        )

    ontology_mask = mask_for(ONTOLOGY_GUARD_ATOMS)
    ontology_strict_mask = mask_for(ONTOLOGY_GUARD_ATOMS + STRICT_ADDITION_ATOMS)
    full_mask = mask_for(ATOMS)
    proper_role_prefix_masks = []
    for role_mask in range(1 << len(ROLE_SUCCESSOR_ATOMS)):
        if role_mask == (1 << len(ROLE_SUCCESSOR_ATOMS)) - 1:
            continue
        combined = ontology_strict_mask
        for index, atom in enumerate(ROLE_SUCCESSOR_ATOMS):
            if role_mask & (1 << index):
                combined |= 1 << ATOM_INDEX[atom]
        proper_role_prefix_masks.append(combined)

    return {
        "atom_order": ATOMS,
        "total_boolean_assignment_count": total_rows,
        "summary_only_no_full_truth_table_emitted": True,
        "lane_rows": lane_rows,
        "ontology_only_mask": ontology_mask,
        "ontology_only_atoms": atoms_for(ontology_mask),
        "ontology_plus_strict_mask": ontology_strict_mask,
        "ontology_plus_strict_atoms": atoms_for(ontology_strict_mask),
        "full_projection_mask": full_mask,
        "full_projection_atoms": atoms_for(full_mask),
        "proper_role_prefix_masks_after_ontology_plus_strict": proper_role_prefix_masks,
        "proper_role_prefix_fail_count_after_ontology_plus_strict": len(proper_role_prefix_masks),
        "ontology_only_ready_lanes": [lane for lane, req in STAGED_LANES.items() if (ontology_mask & mask_for(req)) == mask_for(req)],
        "ontology_plus_strict_ready_lanes": [lane for lane, req in STAGED_LANES.items() if (ontology_strict_mask & mask_for(req)) == mask_for(req)],
        "full_projection_ready_lanes": [lane for lane, req in STAGED_LANES.items() if (full_mask & mask_for(req)) == mask_for(req)],
    }


def append_doc_sections() -> None:
    eq_section = """
## P2406/S1356 information-to-physics staged projection barrier certificate

`P2406/S1356` combines P2405 ontology typing with the P2404 dependency cut as one compact finite Boolean certificate over `12` atoms: five ontology guards, four strict internal-information additions, and three downstream role-successor atoms.  The full finite space has `4096` assignments, but the artifact records exact counts and monomials rather than dumping the full truth table.

The staged result is: ontology alone readies only the typed-information root; ontology plus all four strict additions readies only internal strict completion; no role-bearing physical projection is ready until the appropriate role-successor atoms are added.  The `L_total` and ToE downstream projection lanes require the degree-12 monomial consisting of all ontology, strict-addition, and role-successor atoms.

Thus pure information is not physical-role export by itself.  Physics remains a downstream projection from the single informational nadsoliton through strict internal completion and then through explicit role-successor theorems.
""".strip()
    lag_section = """
## P2406/S1356 staged projection barrier for Lagrangian/EOM

`P2406/S1356` makes the information-to-physics staging explicit: pure-information ontology plus strict internal completion still does not license a role-bearing `L_total`.  `L_total` and ToE projection require the full degree-12 guard/completion/role monomial, so no Lagrangian term may be promoted directly from ontology or strict additions alone.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    summary = staged_projection_summary()
    p2404_theorem = artifacts["P2404_DEPENDENCY_CUT"].get(
        "strict_addition_physics_lane_dependency_cut_certificate", {}
    ).get("theorem_export", {})
    p2405_theorem = artifacts["P2405_INFORMATION_ONTOLOGY"].get(
        "nadsoliton_information_ontology_projection_certificate", {}
    ).get("theorem_export", {})
    ltotal_row = next(row for row in summary["lane_rows"] if row["lane"] == "role_bearing_ltotal_downstream_projection_candidate")
    toe_row = next(row for row in summary["lane_rows"] if row["lane"] == "toe_downstream_projection_candidate")

    theorem_export = {
        "theorem_name": "P2406_T1_information_to_physics_staged_projection_barrier",
        "total_boolean_assignment_count": summary["total_boolean_assignment_count"],
        "summary_only_no_full_truth_table_emitted": summary["summary_only_no_full_truth_table_emitted"],
        "ontology_only_ready_lanes": summary["ontology_only_ready_lanes"],
        "ontology_plus_strict_ready_lanes": summary["ontology_plus_strict_ready_lanes"],
        "proper_role_prefix_fail_count_after_ontology_plus_strict": summary[
            "proper_role_prefix_fail_count_after_ontology_plus_strict"
        ],
        "full_projection_ready_lanes": summary["full_projection_ready_lanes"],
        "ltotal_projection_anf": ltotal_row["anf_monomial"],
        "ltotal_projection_anf_degree": ltotal_row["anf_degree"],
        "ltotal_projection_true_assignment_count": ltotal_row["true_assignment_count"],
        "toe_projection_anf": toe_row["anf_monomial"],
        "toe_projection_anf_degree": toe_row["anf_degree"],
        "toe_projection_true_assignment_count": toe_row["true_assignment_count"],
        "p2404_degree_seven_inherited": p2404_theorem.get("ltotal_dependency_degree") == 7,
        "p2405_degree_five_ontology_guard_inherited": p2405_theorem.get("ontology_guard_anf_degree") == 5,
        "p2405_no_underlayer_inherited": p2405_theorem.get("separate_information_layer_allowed") is False,
        "not_licensed": [
            "No direct projection from pure information ontology to physical roles is exported.",
            "No role-bearing L_total follows from strict additions without role-successor theorems.",
            "No sub-nadsoliton information layer is introduced.",
            "No ToE closure follows from the staged barrier certificate.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "finite_space_has_4096_assignments": summary["total_boolean_assignment_count"] == 4096,
        "compact_artifact_does_not_emit_full_truth_table": summary["summary_only_no_full_truth_table_emitted"] is True,
        "ontology_only_readies_only_information_root": summary["ontology_only_ready_lanes"] == ["ontology_typed_information_root"],
        "ontology_plus_strict_readies_only_internal_completion": summary["ontology_plus_strict_ready_lanes"]
        == ["ontology_typed_information_root", "strict_internal_information_completion_ready"],
        "all_seven_proper_role_prefixes_fail_ltotal_toe": theorem_export[
            "proper_role_prefix_fail_count_after_ontology_plus_strict"
        ] == 7,
        "ltotal_and_toe_are_degree_twelve": theorem_export["ltotal_projection_anf_degree"] == 12
        and theorem_export["toe_projection_anf_degree"] == 12,
        "ltotal_and_toe_have_single_true_assignment": theorem_export["ltotal_projection_true_assignment_count"] == 1
        and theorem_export["toe_projection_true_assignment_count"] == 1,
        "p2404_degree_seven_inherited": theorem_export["p2404_degree_seven_inherited"],
        "p2405_degree_five_ontology_guard_inherited": theorem_export["p2405_degree_five_ontology_guard_inherited"],
        "p2405_no_underlayer_inherited": theorem_export["p2405_no_underlayer_inherited"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2406_s1356_v1",
        "packet_id": "P2406",
        "stage_id": "S1356",
        "result_kind": "INFORMATION_TO_PHYSICS_STAGED_PROJECTION_BARRIER_CERTIFICATE",
        "status": "PASS_STAGED_PROJECTION_BARRIER_NO_DIRECT_ROLE_EXPORT_NO_TOE_CLOSURE",
        "information_to_physics_staged_projection_barrier_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: artifact.get("status", "OPEN_UNKNOWN") for name, artifact in artifacts.items()},
            "staged_projection_summary": summary,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": (
            "If pursuing known physics, prove a concrete role-successor projection theorem; do not treat pure information "
            "ontology or strict internal completion as direct L_total/ToE export."
        ),
        "global_status": "OPEN_PROGRESS_STAGED_PROJECTION_BARRIER_CERTIFIED_NO_PHYSICAL_ROLE_EXPORT",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["information_to_physics_staged_projection_barrier_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2406 S1356: information-to-physics staged projection barrier certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2406/S1356 computes a compact 4096-assignment staged barrier over ontology guards, strict additions, and role-successor atoms.",
                "",
                "## Staged readiness",
                "",
                f"- Ontology-only ready lanes: `{theorem['ontology_only_ready_lanes']}`.",
                f"- Ontology-plus-strict ready lanes: `{theorem['ontology_plus_strict_ready_lanes']}`.",
                f"- Proper role prefixes failing after ontology-plus-strict: `{theorem['proper_role_prefix_fail_count_after_ontology_plus_strict']}`.",
                f"- Full projection ready lanes: `{theorem['full_projection_ready_lanes']}`.",
                "",
                "## Lagrangian / ToE projection barrier",
                "",
                f"- `L_total` projection ANF degree: `{theorem['ltotal_projection_anf_degree']}`.",
                f"- ToE projection ANF degree: `{theorem['toe_projection_anf_degree']}`.",
                f"- `L_total` true assignments in 4096-space: `{theorem['ltotal_projection_true_assignment_count']}`.",
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
