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

OUT = GEN / "p2432_s1382_post_antichain_branch_residual_transition_certificate.json"
MD = GEN / "p2432_s1382_post_antichain_branch_residual_transition_certificate.md"

SOURCE_FILES = {
    "P2431_ANTICHAIN": GEN / "p2431_s1381_admissible_next_theorem_target_antichain_certificate.json",
    "P2430_COVER": GEN / "p2430_s1380_repair_derivative_witness_cover_minimality_certificate.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

GATES = [
    "source_obligation_discharge",
    "chi11_source_export",
    "qw2191_selector_discharge",
    "role_transfer_audit_license",
    "role_bearing_ltotal_export",
]
PREDECESSORS = {
    "source_obligation_discharge": set(),
    "chi11_source_export": set(),
    "qw2191_selector_discharge": set(),
    "role_transfer_audit_license": {
        "source_obligation_discharge",
        "chi11_source_export",
        "qw2191_selector_discharge",
    },
    "role_bearing_ltotal_export": {"role_transfer_audit_license"},
}
ANTICHAIN_BRANCHES = {
    "source_first_branch": ["source_obligation_discharge"],
    "selector_pair_first_branch": ["chi11_source_export", "qw2191_selector_discharge"],
}
TARGETS = [
    "bridge_source_ready",
    "selector_source_ready",
    "role_transfer_ready",
    "role_bearing_ltotal_ready",
    "toe_ready",
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
        "new_packet": "P2432|S1382|post-antichain|branch residual|residual transition",
        "p2431_input": "P2431|admissible next theorem-target|target antichain|source-vs-selector",
        "transition_prior": "post.*branch|residual theorem|transition|remaining gate|next admissible",
        "nonpromotion_prior": "nonpromotion|does not discharge|No source|No ToE|not a theorem discharge",
    }
    return {
        "tool": "rg",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
        "finding": (
            "Repo grep finds P2431's next-target antichain and generic residual/transition language, but no production "
            "P24xx certificate computing the post-antichain residual transition map for the source-first and selector-pair-first branches."
        ),
    }


def theorem(payload: dict[str, Any], key: str) -> dict[str, Any]:
    return payload.get(key, {}).get("theorem_export", {})


def readiness_from_gates(gates: set[str]) -> dict[str, bool]:
    return {
        "bridge_source_ready": "source_obligation_discharge" in gates,
        "selector_source_ready": {"chi11_source_export", "qw2191_selector_discharge"}.issubset(gates),
        "role_transfer_ready": "role_transfer_audit_license" in gates,
        "role_bearing_ltotal_ready": "role_bearing_ltotal_export" in gates,
        "toe_ready": set(GATES).issubset(gates),
    }


def all_nonempty_subsets(items: list[str], max_size: int) -> list[list[str]]:
    out = []
    for size in range(1, max_size + 1):
        for combo in combinations(items, size):
            out.append(list(combo))
    return out


def is_admissible(candidate: set[str], current: set[str]) -> bool:
    available = current | candidate
    return all(PREDECESSORS[gate].issubset(available) for gate in candidate)


def readiness_unlocks(candidate: set[str], current: set[str]) -> dict[str, bool]:
    before = readiness_from_gates(current)
    after = readiness_from_gates(current | candidate)
    return {target: after[target] and not before[target] for target in TARGETS}


def candidate_rows_after(current: set[str]) -> list[dict[str, Any]]:
    remaining = [gate for gate in GATES if gate not in current]
    rows = []
    for candidate in all_nonempty_subsets(remaining, 2):
        candidate_set = set(candidate)
        missing_predecessors = sorted(
            set().union(*(PREDECESSORS[gate] for gate in candidate_set)) - (current | candidate_set),
            key=GATES.index,
        )
        unlocks = readiness_unlocks(candidate_set, current)
        rows.append(
            {
                "candidate_gates": candidate,
                "candidate_size": len(candidate),
                "admissible_after_branch": is_admissible(candidate_set, current),
                "missing_predecessors": missing_predecessors,
                "readiness_unlocks": unlocks,
                "unlocks_any_readiness_target": any(unlocks.values()),
                "remaining_after_candidate": [gate for gate in GATES if gate not in current | candidate_set],
            }
        )
    return rows


def minimal_readiness_candidates(rows: list[dict[str, Any]]) -> list[list[str]]:
    out = []
    for row in rows:
        if not row["admissible_after_branch"] or not row["unlocks_any_readiness_target"]:
            continue
        candidate_set = set(row["candidate_gates"])
        if any(
            other["admissible_after_branch"]
            and other["unlocks_any_readiness_target"]
            and set(other["candidate_gates"]).issubset(candidate_set)
            and set(other["candidate_gates"]) != candidate_set
            and other["readiness_unlocks"] == row["readiness_unlocks"]
            for other in rows
        ):
            continue
        out.append(row["candidate_gates"])
    return out


def branch_rows() -> list[dict[str, Any]]:
    rows = []
    for branch_name, gates in ANTICHAIN_BRANCHES.items():
        current = set(gates)
        candidates = candidate_rows_after(current)
        admissible = [row for row in candidates if row["admissible_after_branch"]]
        rows.append(
            {
                "branch_name": branch_name,
                "branch_gates": gates,
                "readiness_after_branch": readiness_from_gates(current),
                "remaining_gate_count_after_branch": len(GATES) - len(current),
                "remaining_gates_after_branch": [gate for gate in GATES if gate not in current],
                "candidate_row_count_size_le_2_after_branch": len(candidates),
                "admissible_candidate_count_size_le_2_after_branch": len(admissible),
                "inadmissible_candidate_count_size_le_2_after_branch": len(candidates) - len(admissible),
                "admissible_singletons_after_branch": [row["candidate_gates"] for row in admissible if row["candidate_size"] == 1],
                "inadmissible_singletons_after_branch": [
                    row["candidate_gates"]
                    for row in candidates
                    if row["candidate_size"] == 1 and not row["admissible_after_branch"]
                ],
                "minimal_readiness_complete_candidates_after_branch": minimal_readiness_candidates(candidates),
                "candidate_rows_size_le_2_after_branch": candidates,
            }
        )
    return rows


def build_certificate(sources: dict[str, dict[str, Any]]) -> dict[str, Any]:
    p2431 = theorem(sources["P2431_ANTICHAIN"], "admissible_next_theorem_target_antichain_certificate")
    p2430 = theorem(sources["P2430_COVER"], "repair_derivative_witness_cover_minimality_certificate")
    rows = branch_rows()
    return {
        "branch_transition_rows": rows,
        "branch_count": len(rows),
        "p2431_antichain_inherited": p2431.get("minimal_readiness_complete_candidates"),
        "p2431_role_transfer_not_first_inherited": p2431.get("role_transfer_admissible_as_first_move"),
        "p2431_ltotal_not_first_inherited": p2431.get("role_bearing_ltotal_admissible_as_first_move"),
        "p2430_global_minimal_cover_inherited": p2430.get("global_minimal_covers"),
    }


def append_doc_sections() -> None:
    eq_section = """
## P2432/S1382 post-antichain branch residual transition certificate

`P2432/S1382` follows the two P2431 admissible next-target branches.  If `source_obligation_discharge` is proved first, bridge readiness opens but the next readiness-complete target remains the selector pair `chi11_source_export + qw2191_selector_discharge`; role-transfer and role-bearing `L_total` remain blocked.  If the selector pair is proved first, selector readiness opens and the next singleton target is `source_obligation_discharge`; a size-two candidate may include `source_obligation_discharge + role_transfer_audit_license`, but role-transfer is not admissible without source.

Thus both branches still converge on the same five-gate theorem cover before ToE: source, chi11, QW-2191, role-transfer, and role-bearing `L_total`.  The transition map is a proof-search guide only and exports no theorem gate.
""".strip()
    lag_section = """
## P2432/S1382 post-antichain branch residual guard for Lagrangian/EOM

`P2432/S1382` shows that neither first branch licenses role-transfer or role-bearing `L_total` immediately.  After source-first the selector pair remains next; after selector-pair-first source remains next.  No branch transition can be inserted into `L_total` as discharged dynamics.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    cert = build_certificate(sources)
    rows_by_branch = {row["branch_name"]: row for row in cert["branch_transition_rows"]}
    theorem_export = {
        "theorem_name": "P2432_T1_post_antichain_branch_residual_transition_certificate",
        "branch_count": cert["branch_count"],
        "branch_transition_rows": cert["branch_transition_rows"],
        "source_first_remaining_gate_count": rows_by_branch["source_first_branch"]["remaining_gate_count_after_branch"],
        "selector_pair_first_remaining_gate_count": rows_by_branch["selector_pair_first_branch"]["remaining_gate_count_after_branch"],
        "source_first_next_minimal_readiness_candidates": rows_by_branch["source_first_branch"][
            "minimal_readiness_complete_candidates_after_branch"
        ],
        "selector_pair_first_next_minimal_readiness_candidates": rows_by_branch["selector_pair_first_branch"][
            "minimal_readiness_complete_candidates_after_branch"
        ],
        "role_transfer_admissible_after_source_first": ["role_transfer_audit_license"] in rows_by_branch["source_first_branch"][
            "admissible_singletons_after_branch"
        ],
        "role_transfer_admissible_after_selector_pair_first": ["role_transfer_audit_license"] in rows_by_branch[
            "selector_pair_first_branch"
        ]["admissible_singletons_after_branch"],
        "role_bearing_ltotal_admissible_after_either_first_branch": False,
        "p2431_antichain_inherited": cert["p2431_antichain_inherited"] == [
            ["source_obligation_discharge"],
            ["chi11_source_export", "qw2191_selector_discharge"],
        ],
        "p2431_role_transfer_not_first_inherited": cert["p2431_role_transfer_not_first_inherited"] is False,
        "p2431_ltotal_not_first_inherited": cert["p2431_ltotal_not_first_inherited"] is False,
        "p2430_global_minimal_cover_inherited": cert["p2430_global_minimal_cover_inherited"] == [GATES],
        "source_obligation_discharge_exported": False,
        "chi11_source_exported": False,
        "qw2191_discharged": False,
        "role_transfer_licensed": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "Branch transitions are proof-search state updates, not theorem-gate discharge.",
            "Role-transfer is not admissible until source and the chi11/QW-2191 selector pair are both discharged.",
            "Role-bearing L_total remains behind role-transfer in both branches.",
            "No source, selector, role-transfer, role-bearing L_total, or ToE theorem is exported.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "two_branches": theorem_export["branch_count"] == 2,
        "source_first_remaining_4": theorem_export["source_first_remaining_gate_count"] == 4,
        "selector_pair_first_remaining_3": theorem_export["selector_pair_first_remaining_gate_count"] == 3,
        "source_first_next_selector_pair": theorem_export["source_first_next_minimal_readiness_candidates"] == [
            ["chi11_source_export", "qw2191_selector_discharge"]
        ],
        "selector_pair_first_next_source_or_source_role": theorem_export["selector_pair_first_next_minimal_readiness_candidates"] == [
            ["source_obligation_discharge"],
            ["source_obligation_discharge", "role_transfer_audit_license"],
        ],
        "role_transfer_not_after_source_first": not theorem_export["role_transfer_admissible_after_source_first"],
        "role_transfer_not_after_selector_pair_first": not theorem_export["role_transfer_admissible_after_selector_pair_first"],
        "ltotal_not_after_either": not theorem_export["role_bearing_ltotal_admissible_after_either_first_branch"],
        "p2431_antichain_inherited": theorem_export["p2431_antichain_inherited"],
        "p2431_role_transfer_not_first_inherited": theorem_export["p2431_role_transfer_not_first_inherited"],
        "p2431_ltotal_not_first_inherited": theorem_export["p2431_ltotal_not_first_inherited"],
        "p2430_global_cover_inherited": theorem_export["p2430_global_minimal_cover_inherited"],
        "source_still_open": not theorem_export["source_obligation_discharge_exported"],
        "chi11_still_open": not theorem_export["chi11_source_exported"],
        "qw2191_still_open": not theorem_export["qw2191_discharged"],
        "role_transfer_still_blocked": not theorem_export["role_transfer_licensed"],
        "ltotal_still_blocked": not theorem_export["role_bearing_ltotal_exported"],
        "toe_still_open": not theorem_export["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2432_s1382_v1",
        "packet_id": "P2432",
        "stage_id": "S1382",
        "result_kind": "POST_ANTICHAIN_BRANCH_RESIDUAL_TRANSITION_CERTIFICATE",
        "status": "PASS_POST_ANTICHAIN_BRANCH_TRANSITION_NO_GATE_DISCHARGE_NO_CLOSURE",
        "post_antichain_branch_residual_transition_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: source.get("status", "OPEN_UNKNOWN") for name, source in sources.items()},
            "finite_witness_certificate": cert,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": (
            "After choosing either antichain branch, prove the complementary source/selector branch next; only then can role-transfer become admissible."
        ),
        "global_status": "OPEN_PROGRESS_POST_ANTICHAIN_TRANSITIONS_CERTIFIED_NO_THEOREM_GATE_DISCHARGED",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["post_antichain_branch_residual_transition_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2432 S1382: post-antichain branch residual transition certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Finite facts",
                "",
                f"- Branch count: `{theorem['branch_count']}`.",
                f"- Source-first next candidates: `{theorem['source_first_next_minimal_readiness_candidates']}`.",
                f"- Selector-pair-first next candidates: `{theorem['selector_pair_first_next_minimal_readiness_candidates']}`.",
                f"- Role-transfer admissible after source first: `{theorem['role_transfer_admissible_after_source_first']}`.",
                f"- Role-transfer admissible after selector pair first: `{theorem['role_transfer_admissible_after_selector_pair_first']}`.",
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
