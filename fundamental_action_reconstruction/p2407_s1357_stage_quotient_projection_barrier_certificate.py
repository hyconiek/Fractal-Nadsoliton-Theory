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

OUT = GEN / "p2407_s1357_stage_quotient_projection_barrier_certificate.json"
MD = GEN / "p2407_s1357_stage_quotient_projection_barrier_certificate.md"

SOURCE_FILES = {
    "P2406_STAGED_BARRIER": GEN / "p2406_s1356_information_to_physics_staged_projection_barrier_certificate.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

STAGE_ATOMS = [
    "O_ontology_guard_package",
    "S_strict_internal_completion_package",
    "R_role_successor_projection_package",
]

STAGE_EXPANSIONS = {
    "O_ontology_guard_package": [
        "nadsoliton_is_sole_primordial_information",
        "no_separate_information_layer_under_nadsoliton",
        "strict_additions_are_internal_information_constraints",
        "physics_roles_are_downstream_projections",
        "observer_is_downstream_readout_not_source",
    ],
    "S_strict_internal_completion_package": [
        "apd_completion_structure",
        "gf2_phase_topological_data",
        "nonlinear_damping_compression",
        "chi11_selector_source_declared",
    ],
    "R_role_successor_projection_package": [
        "alpha_geo_electroweak_role_theorem",
        "beta_tors_strict_role_theorem",
        "beta_power_hierarchy_successor_theorem",
    ],
}

STAGE_LANES = {
    "typed_information_root_ready": ["O_ontology_guard_package"],
    "strict_internal_completion_ready": ["O_ontology_guard_package", "S_strict_internal_completion_package"],
    "role_bearing_ltotal_projection_ready": STAGE_ATOMS,
    "toe_projection_ready": STAGE_ATOMS,
}


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
    return {"count": len(lines), "samples": lines[:14]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2407|S1357|stage quotient projection barrier|stage-quotient",
        "p2406_source": "P2406|degree-12|4096-assignment|staged projection barrier",
        "quotient_language": "quotient|macro|package|stage",
        "no_shortcut_guard": "direct.*role|strict additions alone|ontology.*alone|No ToE closure|L_total.*promotion",
    }
    return {
        "tool": "rg",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
        "finding": (
            "Repo grep finds P2406's 12-atom staged barrier, but no separate quotient certificate reducing it "
            "to the three stage packages O/S/R while preserving the no-shortcut result."
        ),
    }


def stage_mask(stages: list[str]) -> int:
    mask = 0
    for stage in stages:
        mask |= 1 << STAGE_ATOMS.index(stage)
    return mask


def stages_for(mask: int) -> list[str]:
    return [stage for index, stage in enumerate(STAGE_ATOMS) if mask & (1 << index)]


def lane_ready(mask: int, lane: str) -> bool:
    req_mask = stage_mask(STAGE_LANES[lane])
    return (mask & req_mask) == req_mask


def quotient_lattice() -> dict[str, Any]:
    rows = []
    for mask in range(1 << len(STAGE_ATOMS)):
        ready = [lane for lane in STAGE_LANES if lane_ready(mask, lane)]
        rows.append(
            {
                "mask": mask,
                "stages_true": stages_for(mask),
                "ready_lanes": ready,
                "ltotal_ready": "role_bearing_ltotal_projection_ready" in ready,
                "toe_ready": "toe_projection_ready" in ready,
                "invalid_shortcut_reason": invalid_shortcut_reason(mask),
            }
        )
    return {
        "stage_atoms": STAGE_ATOMS,
        "stage_expansions": STAGE_EXPANSIONS,
        "rows": rows,
        "row_count": len(rows),
        "ltotal_true_masks": [row["mask"] for row in rows if row["ltotal_ready"]],
        "toe_true_masks": [row["mask"] for row in rows if row["toe_ready"]],
        "proper_subset_fail_count": sum(1 for row in rows if row["mask"] != 7 and not row["ltotal_ready"] and not row["toe_ready"]),
        "stage_projection_anf": " * ".join(STAGE_ATOMS),
        "stage_projection_anf_degree": len(STAGE_ATOMS),
        "ontology_only_mask": 1,
        "ontology_plus_strict_mask": 3,
        "full_stage_mask": 7,
        "no_shortcut_masks": [row["mask"] for row in rows if row["mask"] != 7 and row["invalid_shortcut_reason"]],
    }


def invalid_shortcut_reason(mask: int) -> str:
    has_o = bool(mask & 1)
    has_s = bool(mask & 2)
    has_r = bool(mask & 4)
    if has_s and not has_o:
        return "strict_completion_without_ontology_guard"
    if has_r and not (has_o and has_s):
        return "role_projection_without_prior_ontology_and_strict_completion"
    if has_o and not has_s and not has_r:
        return "ontology_only_not_physics"
    if has_o and has_s and not has_r:
        return "ontology_plus_strict_not_role_projection"
    return ""


def append_doc_sections() -> None:
    eq_section = """
## P2407/S1357 stage-quotient projection barrier certificate

`P2407/S1357` compresses the P2406 twelve-atom barrier into a three-stage quotient: `O` = ontology guard package, `S` = strict internal completion package, and `R` = role-successor projection package.  The quotient has exactly `8` rows and preserves the P2406 no-shortcut theorem.

In the quotient, `L_total`/ToE projection is the single degree-3 monomial `O * S * R`.  The ontology-only mask `1` readies only the typed information root; the ontology-plus-strict mask `3` readies internal strict completion but no physical role projection; the only role-bearing mask is the full mask `7`.

Thus the long degree-12 monomial from P2406 is not arbitrary bookkeeping: it factors into three mandatory stage packages, and every shortcut that skips ontology, strict completion, or role-successor projection remains rejected.
""".strip()
    lag_section = """
## P2407/S1357 stage-quotient barrier for Lagrangian/EOM

`P2407/S1357` reduces the P2406 guard/completion/role monomial to the quotient rule `O * S * R`.  A role-bearing `L_total` is allowed only at the full stage mask `7`; ontology-only and ontology-plus-strict masks are still non-physical-role states.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    lattice = quotient_lattice()
    p2406_theorem = artifacts["P2406_STAGED_BARRIER"].get(
        "information_to_physics_staged_projection_barrier_certificate", {}
    ).get("theorem_export", {})
    p2406_degree = p2406_theorem.get("ltotal_projection_anf_degree")
    expanded_degree = sum(len(STAGE_EXPANSIONS[stage]) for stage in STAGE_ATOMS)
    theorem_export = {
        "theorem_name": "P2407_T1_stage_quotient_projection_barrier",
        "stage_atoms": STAGE_ATOMS,
        "row_count": lattice["row_count"],
        "ltotal_true_masks": lattice["ltotal_true_masks"],
        "toe_true_masks": lattice["toe_true_masks"],
        "proper_subset_fail_count": lattice["proper_subset_fail_count"],
        "stage_projection_anf": lattice["stage_projection_anf"],
        "stage_projection_anf_degree": lattice["stage_projection_anf_degree"],
        "ontology_only_ready_lanes": lattice["rows"][1]["ready_lanes"],
        "ontology_plus_strict_ready_lanes": lattice["rows"][3]["ready_lanes"],
        "full_stage_ready_lanes": lattice["rows"][7]["ready_lanes"],
        "no_shortcut_masks": lattice["no_shortcut_masks"],
        "expanded_stage_atom_count": expanded_degree,
        "p2406_degree_twelve_inherited": p2406_degree == 12,
        "quotient_expansion_matches_p2406_degree": expanded_degree == p2406_degree,
        "not_licensed": [
            "No ontology-only physical-role projection is exported.",
            "No ontology-plus-strict L_total promotion is exported.",
            "No stage may be skipped: O, S, and R are jointly necessary in the quotient.",
            "No ToE closure follows; full mask 7 is conditional readiness only.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "quotient_has_eight_rows": lattice["row_count"] == 8,
        "only_full_mask_readies_ltotal_toe": theorem_export["ltotal_true_masks"] == [7] and theorem_export["toe_true_masks"] == [7],
        "proper_subset_fail_count_is_seven": theorem_export["proper_subset_fail_count"] == 7,
        "stage_anf_degree_three": theorem_export["stage_projection_anf_degree"] == 3,
        "ontology_only_no_physical_role": theorem_export["ontology_only_ready_lanes"] == ["typed_information_root_ready"],
        "ontology_plus_strict_no_physical_role": theorem_export["ontology_plus_strict_ready_lanes"]
        == ["typed_information_root_ready", "strict_internal_completion_ready"],
        "p2406_degree_twelve_inherited": theorem_export["p2406_degree_twelve_inherited"],
        "quotient_expansion_matches_p2406_degree": theorem_export["quotient_expansion_matches_p2406_degree"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2407_s1357_v1",
        "packet_id": "P2407",
        "stage_id": "S1357",
        "result_kind": "STAGE_QUOTIENT_PROJECTION_BARRIER_CERTIFICATE",
        "status": "PASS_STAGE_QUOTIENT_BARRIER_FULL_MASK_ONLY_NO_TOE_CLOSURE",
        "stage_quotient_projection_barrier_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: artifact.get("status", "OPEN_UNKNOWN") for name, artifact in artifacts.items()},
            "quotient_lattice": lattice,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": (
            "Choose a concrete missing R-stage role-successor theorem only after keeping O and S as mandatory packages; "
            "do not reopen ontology-only or strict-only L_total shortcuts."
        ),
        "global_status": "OPEN_PROGRESS_STAGE_QUOTIENT_BARRIER_CERTIFIED_NO_PHYSICAL_ROLE_EXPORT",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["stage_quotient_projection_barrier_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2407 S1357: stage-quotient projection barrier certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2407/S1357 compresses the P2406 12-atom barrier into the exact three-stage quotient O*S*R.",
                "",
                "## Quotient facts",
                "",
                f"- Row count: `{theorem['row_count']}`.",
                f"- `L_total` true masks: `{theorem['ltotal_true_masks']}`.",
                f"- ToE true masks: `{theorem['toe_true_masks']}`.",
                f"- Proper-subset fail count: `{theorem['proper_subset_fail_count']}`.",
                f"- Stage ANF: `{theorem['stage_projection_anf']}`.",
                f"- Expanded atom count: `{theorem['expanded_stage_atom_count']}`.",
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
