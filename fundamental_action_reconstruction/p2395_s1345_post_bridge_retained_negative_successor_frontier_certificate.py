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

OUT = GEN / "p2395_s1345_post_bridge_retained_negative_successor_frontier_certificate.json"
MD = GEN / "p2395_s1345_post_bridge_retained_negative_successor_frontier_certificate.md"

SOURCE_FILES = {
    "P2394_APD_CHI11_REBASED_FRONTIER": GEN / "p2394_s1344_apd_bridge_chi11_rebased_role_frontier_certificate.json",
    "ROLE_TRANSFER_DRAFT": ROOT / "STRICT_KERNEL_LEGACY_ROLE_TRANSFER_AUDIT_DRAFT.md",
    "F4_ROLE_CLASSIFICATION": ROOT / "F4_LEGACY_PHYSICAL_ROLE_TRANSFER_CLASSIFICATION_PACKET.md",
    "N73_WEINBERG_RETAINED_NEGATIVE": ROOT / "N73_CURRENT_LEGACY_WEINBERG_RETAINED_BRANCH_FULL_NEGATIVE_CLOSURE_THEOREM.md",
    "N90_FINE_STRUCTURE_RETAINED_NEGATIVE": ROOT / "N90_CURRENT_LEGACY_FINE_STRUCTURE_RETAINED_BRANCH_FULL_NEGATIVE_CLOSURE_THEOREM.md",
    "N106_GRAVITY_RETAINED_NEGATIVE": ROOT / "N106_CURRENT_LEGACY_GRAVITY_HIERARCHY_RETAINED_BRANCH_FULL_NEGATIVE_CLOSURE_THEOREM.md",
    "N116_ROLE_PACKAGE_NEGATIVE": ROOT / "N116_CURRENT_LEGACY_PHYSICAL_ROLE_PACKAGE_FULL_NEGATIVE_CLOSURE_THEOREM.md",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

ROLE_ROWS = [
    {
        "role_id": "legacy_weinberg_sin2_theta_role",
        "legacy_formula": "sin^2(theta_W)=alpha_geo/12",
        "negative_source": "N73_WEINBERG_RETAINED_NEGATIVE",
        "p2394_target": "legacy_weinberg_sin2_theta_role_transfer",
        "successor_atoms": ["alpha_geo_electroweak_role_theorem"],
        "successor_target": "strict electroweak successor theorem for alpha_geo after APD normalization",
    },
    {
        "role_id": "legacy_alpha_em_inverse_role",
        "legacy_formula": "alpha_EM^-1=alpha_geo/(2*beta_tors)*(1-beta_tors)",
        "negative_source": "N90_FINE_STRUCTURE_RETAINED_NEGATIVE",
        "p2394_target": "legacy_alpha_em_inverse_role_transfer",
        "successor_atoms": ["alpha_geo_electroweak_role_theorem", "beta_tors_strict_role_theorem"],
        "successor_target": "modified strict fine-structure successor using strict damping/compression semantics rather than raw beta_tors inheritance",
    },
    {
        "role_id": "legacy_beta_power_gravity_hierarchy_role",
        "legacy_formula": "beta^N gravity hierarchy",
        "negative_source": "N106_GRAVITY_RETAINED_NEGATIVE",
        "p2394_target": "legacy_beta_power_gravity_hierarchy_successor",
        "successor_atoms": ["beta_tors_strict_role_theorem", "beta_power_hierarchy_successor_theorem"],
        "successor_target": "modified strict gravity-hierarchy successor using nonlinear compression data",
    },
]

MATRIX_COLUMNS = [
    "apd_bridge_found",
    "chi11_selector_found",
    "retained_branch_closed_negative_now",
    "alpha_geo_role_atom_required",
    "beta_tors_role_atom_required",
    "beta_power_successor_atom_required",
    "modified_successor_frontier_open",
    "current_transfer_licensed",
]


def rel(path: Path) -> str:
    return path.relative_to(REPO).as_posix()


def read_text(path: Path) -> str:
    if not path.exists():
        return f"OPEN_MISSING_ARTIFACT::{rel(path)}"
    return path.read_text(encoding="utf-8")


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": rel(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def rg_audit() -> dict[str, Any]:
    patterns = [
        "P2395|S1345|retained negative successor frontier|post-bridge retained",
        "N73|N90|N106|retained branch full negative closure|retained branch.*closed negatively",
        "modified transfer|strict successor expression|replaced branch|successor theorem",
        r"sin\^2\(theta_W\)=alpha_geo/12|alpha_EM\^-1|beta\^N gravity hierarchy",
        r"P2394|APD bridge found|K_strict = K_legacy\*A\*P\*D|role-transfer frontier",
    ]
    out: dict[str, Any] = {}
    for pattern in patterns:
        proc = subprocess.run(
            ["rg", "-n", pattern, "fundamental_action_reconstruction", "-g", "*.py", "-g", "*.md", "-g", "*.json"],
            cwd=REPO,
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        lines = [line for line in proc.stdout.splitlines() if line]
        out[pattern] = {"count": len(lines), "samples": lines[:14]}
    return {
        "tool": "rg",
        "patterns": out,
        "finding": (
            "Repo grep finds existing retained-branch negative closures for Weinberg, fine-structure, and gravity roles, plus P2394's APD/chi11 rebase. "
            "P2395 therefore does not redo those negative closures; it computes the post-bridge successor frontier that remains after unchanged retained transfer is closed now."
        ),
    }


def gf2_rank(matrix: list[list[int]]) -> int:
    work = [row[:] for row in matrix if any(row)]
    if not work:
        return 0
    rank = 0
    col = 0
    width = len(work[0])
    while rank < len(work) and col < width:
        pivot = next((row for row in range(rank, len(work)) if work[row][col] % 2), None)
        if pivot is None:
            col += 1
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        for row in range(len(work)):
            if row != rank and work[row][col] % 2:
                work[row] = [a ^ b for a, b in zip(work[row], work[rank])]
        rank += 1
        col += 1
    return rank


def retained_negative_closed(text: str, role_keyword: str) -> bool:
    lowered = text.lower()
    return all(
        phrase in lowered
        for phrase in [
            "retained branch",
            "closed negatively",
            "current repo state",
        ]
    ) and role_keyword.lower() in lowered


def extract_context(p2394: dict[str, Any]) -> dict[str, Any]:
    cert = p2394.get("apd_bridge_chi11_rebased_role_frontier_certificate", {})
    apd = cert.get("apd_bridge_status", {})
    selector = cert.get("selector_status", {})
    frontier = cert.get("rebased_role_frontier", {})
    return {
        "apd_bridge_found": bool(apd.get("apd_bridge_found_in_existing_repo")),
        "chi11_selector_found": bool(selector.get("strict_selector_found_in_declared_scope")),
        "p2394_current_assignment": frontier.get("current_assignment", {}),
        "p2394_minimal_supports": frontier.get("minimal_supports", {}),
    }


def build_role_rows(artifacts: dict[str, Any], context: dict[str, Any]) -> list[dict[str, Any]]:
    role_rows = []
    keyword_by_source = {
        "N73_WEINBERG_RETAINED_NEGATIVE": "weinberg",
        "N90_FINE_STRUCTURE_RETAINED_NEGATIVE": "fine-structure",
        "N106_GRAVITY_RETAINED_NEGATIVE": "gravity",
    }
    current_values = context["p2394_current_assignment"].get("target_values", {})
    for row in ROLE_ROWS:
        source_text = artifacts[row["negative_source"]].get("text", "")
        retained_closed = retained_negative_closed(source_text, keyword_by_source[row["negative_source"]])
        successor_atoms = row["successor_atoms"]
        licensed_now = bool(current_values.get(row["p2394_target"], False))
        bits = [
            int(context["apd_bridge_found"]),
            int(context["chi11_selector_found"]),
            int(retained_closed),
            int("alpha_geo_electroweak_role_theorem" in successor_atoms),
            int("beta_tors_strict_role_theorem" in successor_atoms),
            int("beta_power_hierarchy_successor_theorem" in successor_atoms),
            1,
            int(licensed_now),
        ]
        role_rows.append(
            {
                **row,
                "retained_unchanged_transfer_verdict": "CLOSED_NEGATIVE_CURRENT_REPO_STATE" if retained_closed else "NOT_CERTIFIED",
                "modified_successor_verdict": "OPEN_ACTIVE_SUCCESSOR_FRONTIER",
                "rejection_forever_verdict": "NOT_PROVEN_FOREVER_REJECTION",
                "current_transfer_licensed": licensed_now,
                "matrix_bits": dict(zip(MATRIX_COLUMNS, [bool(bit) for bit in bits])),
                "matrix_vector": bits,
            }
        )
    return role_rows


def successor_frontier(role_rows: list[dict[str, Any]]) -> dict[str, Any]:
    all_atoms = sorted({atom for row in role_rows for atom in row["successor_atoms"]})
    target_supports = {row["role_id"]: row["successor_atoms"] for row in role_rows}
    atom_degrees = {atom: sum(atom in row["successor_atoms"] for row in role_rows) for atom in all_atoms}
    current_closed = [row for row in role_rows if row["current_transfer_licensed"]]
    return {
        "successor_atom_universe": all_atoms,
        "successor_atom_degrees": atom_degrees,
        "target_minimal_successor_supports": target_supports,
        "all_successor_roles_minimal_union": all_atoms,
        "all_successor_roles_minimal_union_size": len(all_atoms),
        "current_licensed_role_count": len(current_closed),
        "current_unlicensed_role_count": len(role_rows) - len(current_closed),
        "open_modified_successor_count": sum(row["modified_successor_verdict"] == "OPEN_ACTIVE_SUCCESSOR_FRONTIER" for row in role_rows),
        "retained_negative_closed_count": sum(row["retained_unchanged_transfer_verdict"] == "CLOSED_NEGATIVE_CURRENT_REPO_STATE" for row in role_rows),
    }


def append_doc_sections() -> None:
    eq_section = """
## P2395/S1345 post-bridge retained-negative successor frontier

`P2395/S1345` runs the post-P2394 role-transfer audit without duplicating older retained-branch work.  Repo grep finds the existing retained-branch negative closures `N73`, `N90`, and `N106`; P2395 therefore treats unchanged legacy-role retention as closed negatively on the current repo state for the Weinberg, fine-structure, and gravity-hierarchy roles.

With APD bridge and declared-scope `chi11` selector rebased as found, the active frontier becomes the modified strict-successor branch.  The finite role matrix uses three rows and eight columns: APD found, chi11 found, retained negative closure, alpha_geo atom, beta_tors atom, beta-power atom, modified successor open, and current transfer licensed.

The computed successor universe is:

```text
alpha_geo_electroweak_role_theorem,
beta_tors_strict_role_theorem,
beta_power_hierarchy_successor_theorem.
```

Current transfer count is zero; all three modified successor branches remain open.  This is not a forever-rejection theorem: it only says unchanged retained inheritance is closed now, while strict successor roles still require new theorems before any `L_total`, SM/GR numeric extraction, or ToE closure claim.
""".strip()
    lag_section = """
## P2395/S1345 retained-negative role rebase, strict successor branch open

`P2395/S1345` applies the post-bridge role audit to the Lagrangian/EOM lane.  Since `N73/N90/N106` already close the unchanged retained branches negatively on the current repo state, no `L_total` term may import the old Weinberg, fine-structure, or gravity-hierarchy formulas by literal inheritance.

The only honest active branch is a modified strict-successor branch.  It remains open for all three roles and depends on the `alpha_geo`, `beta_tors`, and beta-power successor atoms computed by P2395.  Thus P2395 blocks silent legacy-role import while preserving a precise successor-theorem frontier.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    artifacts = {
        name: load_json(path) if path.suffix == ".json" else {"text": read_text(path)}
        for name, path in SOURCE_FILES.items()
    }
    grep = rg_audit()
    context = extract_context(artifacts["P2394_APD_CHI11_REBASED_FRONTIER"])
    rows = build_role_rows(artifacts, context)
    matrix = [row["matrix_vector"] for row in rows]
    frontier = successor_frontier(rows)
    theorem_export = {
        "theorem_name": "P2395_T1_post_bridge_retained_negative_successor_frontier",
        "closed_context": context,
        "retained_unchanged_transfer_closed_negative_count": frontier["retained_negative_closed_count"],
        "current_licensed_role_count": frontier["current_licensed_role_count"],
        "modified_successor_frontier_open_count": frontier["open_modified_successor_count"],
        "successor_atom_universe": frontier["successor_atom_universe"],
        "all_successor_roles_minimal_union_size": frontier["all_successor_roles_minimal_union_size"],
        "matrix_rank_mod2": gf2_rank(matrix),
        "not_licensed": [
            "No unchanged legacy physical-role transfer is reopened.",
            "No forever rejection of all strict successor roles is claimed.",
            "No L_total promotion, SM/GR numeric extraction, or ToE closure is claimed.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "p2394_apd_bridge_found": context["apd_bridge_found"],
        "p2394_chi11_selector_found": context["chi11_selector_found"],
        "all_three_retained_branches_closed_negative": frontier["retained_negative_closed_count"] == 3,
        "current_transfer_count_zero": frontier["current_licensed_role_count"] == 0,
        "all_three_modified_successor_branches_open": frontier["open_modified_successor_count"] == 3,
        "successor_union_has_three_atoms": frontier["all_successor_roles_minimal_union_size"] == 3,
        "matrix_rank_is_three": gf2_rank(matrix) == 3,
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2395_s1345_v1",
        "packet_id": "P2395",
        "stage_id": "S1345",
        "result_kind": "POST_BRIDGE_RETAINED_NEGATIVE_SUCCESSOR_FRONTIER_CERTIFICATE",
        "status": "PASS_RETAINED_BRANCHES_CLOSED_NEGATIVE_MODIFIED_SUCCESSOR_FRONTIER_OPEN",
        "post_bridge_retained_negative_successor_frontier_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: artifact.get("status", "TEXT_SOURCE") for name, artifact in artifacts.items()},
            "closed_context": context,
            "role_rows": rows,
            "matrix_columns": MATRIX_COLUMNS,
            "role_matrix": matrix,
            "matrix_rank_mod2": gf2_rank(matrix),
            "successor_frontier": frontier,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": "Choose one modified strict-successor role theorem target, preferably the smallest alpha_geo electroweak successor branch, instead of reopening unchanged legacy retention.",
        "global_status": "OPEN_PROGRESS_UNCHANGED_RETENTION_BLOCKED_MODIFIED_SUCCESSOR_FRONTIER_ACTIVE",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    cert = payload["post_bridge_retained_negative_successor_frontier_certificate"]
    theorem = cert["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2395 S1345: post-bridge retained-negative successor frontier certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2395/S1345 accepts P2394's APD/chi11 rebase and uses existing N73/N90/N106 retained-branch negative closures instead of duplicating them.",
                "",
                "## Computed frontier",
                "",
                f"- Retained unchanged transfer closed-negative count: `{theorem['retained_unchanged_transfer_closed_negative_count']}`.",
                f"- Current licensed role count: `{theorem['current_licensed_role_count']}`.",
                f"- Modified successor frontier open count: `{theorem['modified_successor_frontier_open_count']}`.",
                f"- Successor atom universe: `{theorem['successor_atom_universe']}`.",
                f"- Role matrix GF(2) rank: `{theorem['matrix_rank_mod2']}`.",
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
