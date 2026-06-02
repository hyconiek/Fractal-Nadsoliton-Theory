#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from itertools import combinations
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2434_s1384_conditional_legacy_role_transfer_claim_lattice_certificate.json"
MD = GEN / "p2434_s1384_conditional_legacy_role_transfer_claim_lattice_certificate.md"

SOURCE_FILES = {
    "P2433_CONVERGENCE": GEN / "p2433_s1383_source_selector_convergence_role_transfer_gate_certificate.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

ROLE_OBLIGATIONS = [
    "role_transfer_audit_license",
    "role_bearing_ltotal_export",
    "alpha_geo_strict_role_successor_theorem",
    "beta_tors_strict_role_successor_theorem",
    "strict_nonlinear_hierarchy_successor_theorem",
    "chi11_orientation_role_successor_theorem",
]

LEGACY_ROLE_CLAIMS = [
    {
        "role_id": "legacy_weak_mixing_angle",
        "legacy_formula": "sin^2(theta_W)=alpha_geo/12",
        "required_obligations": [
            "role_transfer_audit_license",
            "role_bearing_ltotal_export",
            "alpha_geo_strict_role_successor_theorem",
        ],
        "blocked_reason": "alpha_geo is a bridge normalization datum; no strict alpha_geo electroweak-role successor theorem is exported.",
    },
    {
        "role_id": "legacy_inverse_alpha_em",
        "legacy_formula": "alpha_EM^-1=alpha_geo/(2*beta_tors)*(1-beta_tors)",
        "required_obligations": [
            "role_transfer_audit_license",
            "role_bearing_ltotal_export",
            "alpha_geo_strict_role_successor_theorem",
            "beta_tors_strict_role_successor_theorem",
        ],
        "blocked_reason": "the legacy formula jointly uses alpha_geo and beta_tors, but the strict bridge exports neither unchanged physical role.",
    },
    {
        "role_id": "legacy_beta_power_gravity_hierarchy",
        "legacy_formula": "beta^N gravity hierarchy / beta_tors-power hierarchy",
        "required_obligations": [
            "role_transfer_audit_license",
            "role_bearing_ltotal_export",
            "beta_tors_strict_role_successor_theorem",
            "strict_nonlinear_hierarchy_successor_theorem",
        ],
        "blocked_reason": "strict nonlinear damping/compression is not a literal beta_tors-power hierarchy theorem.",
    },
    {
        "role_id": "legacy_torsion_to_chi11_orientation",
        "legacy_formula": "candidate beta_tors -> chi_11 orientation/torsion role",
        "required_obligations": [
            "role_transfer_audit_license",
            "role_bearing_ltotal_export",
            "beta_tors_strict_role_successor_theorem",
            "chi11_orientation_role_successor_theorem",
        ],
        "blocked_reason": "P2433 can make chi11/QW-2191 selector readiness conditional, but it still does not prove beta_tors -> chi_11 as a role theorem.",
    },
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
    return {"count": len(lines), "samples": lines[:20]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2434|S1384|conditional legacy role-transfer|role transfer claim lattice|legacy role-transfer claim lattice",
        "p2433_input": "P2433|source-selector convergence|role-transfer gate|role_transfer_audit_license",
        "legacy_role_claims": "sin\\^2\\(theta_W\\)|alpha_EM\\^-1|beta\\^N|beta_tors.*chi_11|legacy physical-role",
        "prior_role_audits": "role-transfer pre-audit|claim-specific strict-side partition|legacy role transfer|GF\\(2\\) dependency matrix",
    }
    return {
        "tool": "rg",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
        "finding": (
            "Repo grep finds scratch/current role-transfer pre-audits and claim-specific partition obstructions, but no production P24xx "
            "certificate enumerating the post-P2433 conditional role-claim obligation lattice and proving role-transfer is not automatic."
        ),
    }


def theorem(payload: dict[str, Any], key: str) -> dict[str, Any]:
    return payload.get(key, {}).get("theorem_export", {})


def subset_rows(items: list[str]) -> list[dict[str, Any]]:
    rows = []
    for mask in range(1 << len(items)):
        true_obligations = [item for index, item in enumerate(items) if (mask >> index) & 1]
        true_set = set(true_obligations)
        ready_roles = [
            claim["role_id"] for claim in LEGACY_ROLE_CLAIMS if set(claim["required_obligations"]).issubset(true_set)
        ]
        rows.append(
            {
                "mask": mask,
                "true_obligations": true_obligations,
                "ready_role_claims": ready_roles,
                "ready_role_claim_count": len(ready_roles),
                "all_legacy_roles_ready": len(ready_roles) == len(LEGACY_ROLE_CLAIMS),
            }
        )
    return rows


def minimal_masks_for_claim(rows: list[dict[str, Any]], role_id: str) -> list[list[str]]:
    candidates = [row["true_obligations"] for row in rows if role_id in row["ready_role_claims"]]
    minimal = []
    for candidate in candidates:
        cset = set(candidate)
        if not any(set(other) < cset for other in candidates):
            minimal.append(candidate)
    return sorted(minimal, key=lambda item: (len(item), item))


def role_claim_rows(lattice_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for claim in LEGACY_ROLE_CLAIMS:
        minimal = minimal_masks_for_claim(lattice_rows, claim["role_id"])
        rows.append(
            {
                **claim,
                "minimal_obligation_masks": minimal,
                "minimal_obligation_size": len(minimal[0]) if minimal else None,
                "current_transfer_verdict": "blocked_not_transferred",
                "transferred_by_this_certificate": False,
            }
        )
    return rows


def count_by_ready_role_count(rows: list[dict[str, Any]]) -> dict[str, int]:
    counts = {str(index): 0 for index in range(len(LEGACY_ROLE_CLAIMS) + 1)}
    for row in rows:
        counts[str(row["ready_role_claim_count"])] += 1
    return counts


def append_doc_sections() -> None:
    eq_section = """
## P2434/S1384 conditional legacy role-transfer claim lattice certificate

`P2434/S1384` takes the P2433 source+selector convergence as a hypothetical input and enumerates the remaining role-transfer obligation lattice for the four legacy physical-role claims: Weinberg angle, inverse alpha_EM, beta-power gravity hierarchy, and beta_tors -> chi_11 orientation.  The lattice has 64 masks over six post-convergence obligations; the current mask transfers zero claims.

The role-transfer audit license and role-bearing `L_total` export are necessary but not sufficient: each legacy role still needs its claim-specific strict successor theorem.  Therefore even post-P2433 convergence cannot silently transfer `sin^2(theta_W)=alpha_geo/12`, `alpha_EM^-1`, `beta^N`, or `beta_tors -> chi_11`; the certificate is a conditional audit map, not a role theorem.
""".strip()
    lag_section = """
## P2434/S1384 conditional legacy-role lattice guard for Lagrangian/EOM

`P2434/S1384` blocks the tempting shortcut from role-transfer admissibility to role-bearing terms.  Even after hypothetical source+selector convergence, `L_total` may not contain legacy Weinberg, alpha_EM, gravity-hierarchy, or beta_tors -> chi_11 role terms until the role-transfer audit license, role-bearing `L_total` export, and the relevant claim-specific successor theorem are all actually exported.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    p2433 = theorem(sources["P2433_CONVERGENCE"], "source_selector_convergence_role_transfer_gate_certificate")
    lattice_rows = subset_rows(ROLE_OBLIGATIONS)
    claim_rows = role_claim_rows(lattice_rows)
    all_roles_masks = [row["mask"] for row in lattice_rows if row["all_legacy_roles_ready"]]
    current_row = lattice_rows[0]
    theorem_export = {
        "theorem_name": "P2434_T1_conditional_legacy_role_transfer_claim_lattice_certificate",
        "role_obligation_count": len(ROLE_OBLIGATIONS),
        "role_obligation_names": ROLE_OBLIGATIONS,
        "legacy_role_claim_count": len(LEGACY_ROLE_CLAIMS),
        "role_claim_ids": [claim["role_id"] for claim in LEGACY_ROLE_CLAIMS],
        "role_lattice_assignment_count": len(lattice_rows),
        "ready_role_count_distribution": count_by_ready_role_count(lattice_rows),
        "current_mask": current_row["mask"],
        "current_ready_role_claims": current_row["ready_role_claims"],
        "current_ready_role_claim_count": current_row["ready_role_claim_count"],
        "all_roles_ready_masks": all_roles_masks,
        "all_roles_ready_mask_count": len(all_roles_masks),
        "minimal_masks_by_role_claim": {row["role_id"]: row["minimal_obligation_masks"] for row in claim_rows},
        "minimal_size_by_role_claim": {row["role_id"]: row["minimal_obligation_size"] for row in claim_rows},
        "role_transfer_and_ltotal_are_necessary_for_every_role": all(
            {"role_transfer_audit_license", "role_bearing_ltotal_export"}.issubset(row["required_obligations"])
            for row in claim_rows
        ),
        "role_transfer_and_ltotal_are_not_sufficient_for_any_role": not any(
            set(row["required_obligations"]).issubset({"role_transfer_audit_license", "role_bearing_ltotal_export"})
            for row in claim_rows
        ),
        "p2433_convergence_ready_inherited": p2433.get("readiness_after_convergence") == {
            "bridge_source_ready": True,
            "selector_source_ready": True,
            "role_transfer_ready": False,
            "role_bearing_ltotal_ready": False,
            "toe_ready": False,
        },
        "p2433_role_transfer_next_inherited": p2433.get("role_transfer_admissible_after_source_selector_convergence") is True,
        "source_obligation_discharge_exported_by_this_certificate": False,
        "chi11_source_exported_by_this_certificate": False,
        "qw2191_discharged_by_this_certificate": False,
        "role_transfer_licensed": False,
        "role_bearing_ltotal_exported": False,
        "legacy_role_claim_transferred_by_this_certificate": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "P2434 assumes only a hypothetical post-P2433 convergence state; it does not discharge source or selector gates.",
            "A role-transfer audit license plus role-bearing L_total export is necessary but not sufficient for any legacy role claim.",
            "Every legacy role claim still needs a claim-specific strict successor theorem.",
            "No Weinberg, alpha_EM, beta^N, beta_tors -> chi_11, role-transfer, L_total, or ToE theorem is exported.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "six_obligations": theorem_export["role_obligation_count"] == 6,
        "four_role_claims": theorem_export["legacy_role_claim_count"] == 4,
        "sixty_four_assignments": theorem_export["role_lattice_assignment_count"] == 64,
        "current_transfers_zero": theorem_export["current_ready_role_claim_count"] == 0,
        "one_all_roles_mask": theorem_export["all_roles_ready_mask_count"] == 1 and theorem_export["all_roles_ready_masks"] == [63],
        "weak_mixing_minimal_size_three": theorem_export["minimal_size_by_role_claim"]["legacy_weak_mixing_angle"] == 3,
        "other_roles_minimal_size_four": all(
            theorem_export["minimal_size_by_role_claim"][role_id] == 4
            for role_id in [
                "legacy_inverse_alpha_em",
                "legacy_beta_power_gravity_hierarchy",
                "legacy_torsion_to_chi11_orientation",
            ]
        ),
        "role_transfer_and_ltotal_necessary": theorem_export["role_transfer_and_ltotal_are_necessary_for_every_role"],
        "role_transfer_and_ltotal_not_sufficient": theorem_export["role_transfer_and_ltotal_are_not_sufficient_for_any_role"],
        "p2433_convergence_inherited": theorem_export["p2433_convergence_ready_inherited"],
        "p2433_role_transfer_next_inherited": theorem_export["p2433_role_transfer_next_inherited"],
        "no_source_export": not theorem_export["source_obligation_discharge_exported_by_this_certificate"],
        "no_chi11_export": not theorem_export["chi11_source_exported_by_this_certificate"],
        "no_qw2191_export": not theorem_export["qw2191_discharged_by_this_certificate"],
        "no_role_transfer_export": not theorem_export["role_transfer_licensed"],
        "no_ltotal_export": not theorem_export["role_bearing_ltotal_exported"],
        "no_legacy_role_transfer": not theorem_export["legacy_role_claim_transferred_by_this_certificate"],
        "no_toe_export": not theorem_export["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2434_s1384_v1",
        "packet_id": "P2434",
        "stage_id": "S1384",
        "result_kind": "CONDITIONAL_LEGACY_ROLE_TRANSFER_CLAIM_LATTICE_CERTIFICATE",
        "status": "PASS_CONDITIONAL_LEGACY_ROLE_TRANSFER_CLAIM_LATTICE_NO_ROLE_NO_LTOTAL_NO_TOE_CLOSURE",
        "conditional_legacy_role_transfer_claim_lattice_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: source.get("status", "OPEN_UNKNOWN") for name, source in sources.items()},
            "role_claim_rows": claim_rows,
            "role_lattice_rows": lattice_rows,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": (
            "After real source/selector discharge, prove role-transfer audit license first, then claim-specific strict successor theorems; do not import legacy formulas unchanged."
        ),
        "global_status": "OPEN_PROGRESS_CONDITIONAL_ROLE_TRANSFER_LATTICE_CERTIFIED_NO_LEGACY_ROLE_IMPORTED",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["conditional_legacy_role_transfer_claim_lattice_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2434 S1384: conditional legacy role-transfer claim lattice certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Finite facts",
                "",
                f"- Role obligations: `{theorem['role_obligation_count']}`.",
                f"- Legacy role claims: `{theorem['legacy_role_claim_count']}`.",
                f"- Lattice assignments: `{theorem['role_lattice_assignment_count']}`.",
                f"- Current ready role claims: `{theorem['current_ready_role_claims']}`.",
                f"- Ready-role distribution: `{theorem['ready_role_count_distribution']}`.",
                f"- All-role ready masks: `{theorem['all_roles_ready_masks']}`.",
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
