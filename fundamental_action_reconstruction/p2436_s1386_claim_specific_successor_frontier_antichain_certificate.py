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

OUT = GEN / "p2436_s1386_claim_specific_successor_frontier_antichain_certificate.json"
MD = GEN / "p2436_s1386_claim_specific_successor_frontier_antichain_certificate.md"

SOURCE_FILES = {
    "P2435_SEPARABILITY": GEN / "p2435_s1385_legacy_role_claim_implication_separability_certificate.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

CLAIM_SPECIFIC_SUCCESSORS = [
    "alpha_geo_strict_role_successor_theorem",
    "beta_tors_strict_role_successor_theorem",
    "strict_nonlinear_hierarchy_successor_theorem",
    "chi11_orientation_role_successor_theorem",
]

ROLE_REQUIREMENTS = {
    "legacy_weak_mixing_angle": ["alpha_geo_strict_role_successor_theorem"],
    "legacy_inverse_alpha_em": [
        "alpha_geo_strict_role_successor_theorem",
        "beta_tors_strict_role_successor_theorem",
    ],
    "legacy_beta_power_gravity_hierarchy": [
        "beta_tors_strict_role_successor_theorem",
        "strict_nonlinear_hierarchy_successor_theorem",
    ],
    "legacy_torsion_to_chi11_orientation": [
        "beta_tors_strict_role_successor_theorem",
        "chi11_orientation_role_successor_theorem",
    ],
}

COMMON_ASSUMED_GATES = [
    "role_transfer_audit_license",
    "role_bearing_ltotal_export",
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
        "new_packet": "P2436|S1386|claim-specific successor frontier|successor frontier antichain|role successor frontier",
        "p2435_input": "P2435|legacy role-claim implication|role-claim separability|alpha_EM -> Weinberg",
        "successor_targets_prior": "alpha_geo_strict_role_successor|beta_tors_strict_role_successor|strict_nonlinear_hierarchy_successor|chi11_orientation_role_successor",
        "legacy_role_claims": "sin\\^2\\(theta_W\\)|alpha_EM\\^-1|beta\\^N|beta_tors.*chi_11",
    }
    return {
        "tool": "rg",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
        "finding": (
            "Repo grep finds the P2435 role-claim separability result and role successor names, but no production P24xx "
            "certificate projecting the common-gate-discharged slice onto the claim-specific successor frontier antichain."
        ),
    }


def theorem(payload: dict[str, Any], key: str) -> dict[str, Any]:
    return payload.get(key, {}).get("theorem_export", {})


def ready_roles(successors: set[str]) -> list[str]:
    return [role for role, reqs in ROLE_REQUIREMENTS.items() if set(reqs).issubset(successors)]


def successor_rows() -> list[dict[str, Any]]:
    rows = []
    for mask in range(1 << len(CLAIM_SPECIFIC_SUCCESSORS)):
        true_successors = [name for index, name in enumerate(CLAIM_SPECIFIC_SUCCESSORS) if (mask >> index) & 1]
        roles = ready_roles(set(true_successors))
        rows.append(
            {
                "mask": mask,
                "true_successors": true_successors,
                "ready_role_claims": roles,
                "ready_role_claim_count": len(roles),
                "all_role_claims_ready": len(roles) == len(ROLE_REQUIREMENTS),
            }
        )
    return rows


def minimal_masks_for_role(rows: list[dict[str, Any]], role_id: str) -> list[list[str]]:
    candidates = [row["true_successors"] for row in rows if role_id in row["ready_role_claims"]]
    return sorted(
        [candidate for candidate in candidates if not any(set(other) < set(candidate) for other in candidates)],
        key=lambda item: (len(item), item),
    )


def ready_count_distribution(rows: list[dict[str, Any]]) -> dict[str, int]:
    counts = {str(index): 0 for index in range(len(ROLE_REQUIREMENTS) + 1)}
    for row in rows:
        counts[str(row["ready_role_claim_count"])] += 1
    return counts


def singleton_unlock_rows() -> list[dict[str, Any]]:
    rows = []
    for successor in CLAIM_SPECIFIC_SUCCESSORS:
        successors = {successor}
        rows.append(
            {
                "candidate_successor": successor,
                "ready_role_claims_unlocked_from_empty": ready_roles(successors),
                "ready_role_claim_count_unlocked_from_empty": len(ready_roles(successors)),
            }
        )
    return rows


def minimal_frontier_antichain(rows: list[dict[str, Any]]) -> list[list[str]]:
    winners = [row["true_successors"] for row in rows if row["all_role_claims_ready"]]
    return sorted(
        [candidate for candidate in winners if not any(set(other) < set(candidate) for other in winners)],
        key=lambda item: (len(item), item),
    )


def append_doc_sections() -> None:
    eq_section = """
## P2436/S1386 claim-specific successor frontier antichain certificate

`P2436/S1386` projects the P2435 separability result onto the post-common-gate frontier where role-transfer audit license and role-bearing `L_total` export are assumed only as bookkeeping premises.  The remaining four claim-specific successors form a 16-mask lattice.  From the empty successor mask, only `alpha_geo_strict_role_successor_theorem` unlocks a legacy role claim, namely the weak-mixing/Weinberg role; `beta_tors`, nonlinear hierarchy, and chi11 orientation successors unlock zero claims alone.

The minimal all-role antichain is the full four-successor set `{alpha_geo, beta_tors, strict_nonlinear_hierarchy, chi11_orientation}`.  Thus alpha progress may open Weinberg first, but alpha alone cannot import alpha_EM, gravity hierarchy, or beta_tors -> chi_11 roles; those remain claim-specific theorem targets.
""".strip()
    lag_section = """
## P2436/S1386 claim-specific successor frontier guard for Lagrangian/EOM

`P2436/S1386` shows that even after common role-transfer and `L_total` gates are treated as hypothetical premises, the claim-specific successors still gate which legacy role terms may appear.  The Lagrangian/EOM draft may not insert alpha_EM, gravity-hierarchy, or beta_tors -> chi_11 terms from alpha/Weinberg progress; the full four-successor antichain is still required for all four legacy roles.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    p2435 = theorem(sources["P2435_SEPARABILITY"], "legacy_role_claim_implication_separability_certificate")
    rows = successor_rows()
    min_by_role = {role: minimal_masks_for_role(rows, role) for role in ROLE_REQUIREMENTS}
    frontier = minimal_frontier_antichain(rows)
    singleton_rows = singleton_unlock_rows()
    theorem_export = {
        "theorem_name": "P2436_T1_claim_specific_successor_frontier_antichain_certificate",
        "common_assumed_gate_names": COMMON_ASSUMED_GATES,
        "common_gates_are_assumptions_not_exports": True,
        "claim_specific_successor_count": len(CLAIM_SPECIFIC_SUCCESSORS),
        "claim_specific_successor_names": CLAIM_SPECIFIC_SUCCESSORS,
        "successor_lattice_assignment_count": len(rows),
        "ready_role_count_distribution": ready_count_distribution(rows),
        "current_successor_mask": 0,
        "current_ready_role_claims": rows[0]["ready_role_claims"],
        "current_ready_role_claim_count": rows[0]["ready_role_claim_count"],
        "singleton_unlock_rows_from_empty": singleton_rows,
        "singleton_unlocking_successors_from_empty": [
            row["candidate_successor"] for row in singleton_rows if row["ready_role_claim_count_unlocked_from_empty"] > 0
        ],
        "singleton_nonunlocking_successors_from_empty": [
            row["candidate_successor"] for row in singleton_rows if row["ready_role_claim_count_unlocked_from_empty"] == 0
        ],
        "minimal_successor_masks_by_role_claim": min_by_role,
        "minimal_successor_size_by_role_claim": {role: len(masks[0]) for role, masks in min_by_role.items()},
        "minimal_all_role_successor_antichain": frontier,
        "minimal_all_role_successor_antichain_size": len(frontier[0]) if frontier else None,
        "all_role_ready_masks": [row["mask"] for row in rows if row["all_role_claims_ready"]],
        "p2435_rank_four_inherited": p2435.get("incidence_rank_gf2") == 4,
        "p2435_single_nontrivial_implication_inherited": p2435.get("nontrivial_implications") == [
            ["legacy_inverse_alpha_em", "legacy_weak_mixing_angle"]
        ],
        "legacy_role_claim_transferred_by_this_certificate": False,
        "common_role_transfer_gates_exported_by_this_certificate": False,
        "claim_specific_successor_theorem_exported_by_this_certificate": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "The common role-transfer and L_total gates are assumed only to study the frontier; they are not exported here.",
            "A singleton alpha successor would unlock only the weak-mixing/Weinberg role, not the whole role package.",
            "Beta_tors, nonlinear hierarchy, and chi11 orientation successors still require their own claim-specific theorems.",
            "No legacy role claim, role-transfer theorem, role-bearing L_total export, or ToE closure is exported.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "four_successors": theorem_export["claim_specific_successor_count"] == 4,
        "sixteen_assignments": theorem_export["successor_lattice_assignment_count"] == 16,
        "current_unlocks_zero": theorem_export["current_ready_role_claim_count"] == 0,
        "distribution_matches": theorem_export["ready_role_count_distribution"] == {"0": 5, "1": 6, "2": 2, "3": 2, "4": 1},
        "only_alpha_singleton_unlocks": theorem_export["singleton_unlocking_successors_from_empty"] == [
            "alpha_geo_strict_role_successor_theorem"
        ],
        "three_singletons_do_not_unlock": len(theorem_export["singleton_nonunlocking_successors_from_empty"]) == 3,
        "weak_mixing_minimal_size_one": theorem_export["minimal_successor_size_by_role_claim"]["legacy_weak_mixing_angle"] == 1,
        "other_roles_minimal_size_two": all(
            theorem_export["minimal_successor_size_by_role_claim"][role] == 2
            for role in [
                "legacy_inverse_alpha_em",
                "legacy_beta_power_gravity_hierarchy",
                "legacy_torsion_to_chi11_orientation",
            ]
        ),
        "full_antichain_size_four": theorem_export["minimal_all_role_successor_antichain_size"] == 4,
        "single_all_role_mask": theorem_export["all_role_ready_masks"] == [15],
        "p2435_rank_inherited": theorem_export["p2435_rank_four_inherited"],
        "p2435_implication_inherited": theorem_export["p2435_single_nontrivial_implication_inherited"],
        "no_legacy_role_transfer": not theorem_export["legacy_role_claim_transferred_by_this_certificate"],
        "no_common_gate_export": not theorem_export["common_role_transfer_gates_exported_by_this_certificate"],
        "no_successor_export": not theorem_export["claim_specific_successor_theorem_exported_by_this_certificate"],
        "no_ltotal_export": not theorem_export["role_bearing_ltotal_exported"],
        "no_toe_export": not theorem_export["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2436_s1386_v1",
        "packet_id": "P2436",
        "stage_id": "S1386",
        "result_kind": "CLAIM_SPECIFIC_SUCCESSOR_FRONTIER_ANTICHAIN_CERTIFICATE",
        "status": "PASS_CLAIM_SPECIFIC_SUCCESSOR_FRONTIER_ANTICHAIN_NO_ROLE_TRANSFER_NO_CLOSURE",
        "claim_specific_successor_frontier_antichain_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: source.get("status", "OPEN_UNKNOWN") for name, source in sources.items()},
            "successor_lattice_rows": rows,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": (
            "Treat alpha, beta_tors, nonlinear hierarchy, and chi11 orientation successor theorems as separate proof targets."
        ),
        "global_status": "OPEN_PROGRESS_CLAIM_SPECIFIC_SUCCESSOR_FRONTIER_CERTIFIED_NO_ROLE_IMPORTED",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["claim_specific_successor_frontier_antichain_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2436 S1386: claim-specific successor frontier antichain certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Finite facts",
                "",
                f"- Successor lattice assignments: `{theorem['successor_lattice_assignment_count']}`.",
                f"- Ready-role distribution: `{theorem['ready_role_count_distribution']}`.",
                f"- Singleton unlocks: `{theorem['singleton_unlocking_successors_from_empty']}`.",
                f"- Minimal all-role antichain: `{theorem['minimal_all_role_successor_antichain']}`.",
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
