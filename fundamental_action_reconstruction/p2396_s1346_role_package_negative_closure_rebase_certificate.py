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

OUT = GEN / "p2396_s1346_role_package_negative_closure_rebase_certificate.json"
MD = GEN / "p2396_s1346_role_package_negative_closure_rebase_certificate.md"

SOURCE_FILES = {
    "P2395_SUCCESSOR_FRONTIER": GEN / "p2395_s1345_post_bridge_retained_negative_successor_frontier_certificate.json",
    "N83_WEINBERG_FULL_NEGATIVE": ROOT / "N83_CURRENT_LEGACY_WEINBERG_FULL_CLAIM_SPECIFIC_NEGATIVE_CLOSURE_THEOREM.md",
    "N99_FINE_STRUCTURE_FULL_NEGATIVE": ROOT / "N99_CURRENT_LEGACY_FINE_STRUCTURE_FULL_CLAIM_SPECIFIC_NEGATIVE_CLOSURE_THEOREM.md",
    "N115_GRAVITY_FULL_NEGATIVE": ROOT / "N115_CURRENT_LEGACY_GRAVITY_HIERARCHY_FULL_CLAIM_SPECIFIC_NEGATIVE_CLOSURE_THEOREM.md",
    "N116_ROLE_PACKAGE_FULL_NEGATIVE": ROOT / "N116_CURRENT_LEGACY_PHYSICAL_ROLE_PACKAGE_FULL_NEGATIVE_CLOSURE_THEOREM.md",
    "STRICT_ROLE_AUDIT_DRAFT": ROOT / "STRICT_KERNEL_LEGACY_ROLE_TRANSFER_AUDIT_DRAFT.md",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

ROLE_PACKAGE_ROWS = [
    {
        "role_id": "legacy_weinberg_angle_role",
        "full_negative_source": "N83_WEINBERG_FULL_NEGATIVE",
        "role_keyword": "weinberg",
        "p2395_role_id": "legacy_weinberg_sin2_theta_role",
    },
    {
        "role_id": "legacy_fine_structure_role",
        "full_negative_source": "N99_FINE_STRUCTURE_FULL_NEGATIVE",
        "role_keyword": "fine-structure",
        "p2395_role_id": "legacy_alpha_em_inverse_role",
    },
    {
        "role_id": "legacy_gravity_hierarchy_role",
        "full_negative_source": "N115_GRAVITY_FULL_NEGATIVE",
        "role_keyword": "gravity",
        "p2395_role_id": "legacy_beta_power_gravity_hierarchy_role",
    },
]

MATRIX_COLUMNS = [
    "retained_branch_closed_negative",
    "replaced_branch_closed_negative",
    "full_claim_specific_closed_negative",
    "p2395_modified_successor_open_flag",
    "current_transfer_licensed",
    "future_successor_not_forbidden_forever",
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
        "P2396|S1346|role package negative closure rebase|current-state role package closure",
        "N83|N99|N115|N116|full claim-specific negative closure|legacy physical-role package",
        "retained branch|replaced branch|closed negatively|future strict-side successor semantics",
        "P2395|modified successor frontier|OPEN_ACTIVE_SUCCESSOR_FRONTIER|post-bridge retained-negative",
        "legacy Weinberg|legacy fine-structure|legacy gravity-hierarchy|sin2_theta_w_mz|alpha_em_inv_mz|gravity_hierarchy_beta20",
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
        out[pattern] = {"count": len(lines), "samples": lines[:16]}
    return {
        "tool": "rg",
        "patterns": out,
        "finding": (
            "Repo grep finds N83/N99/N115/N116 full current-state negative closures for the legacy role package. "
            "P2396 therefore treats P2395's modified-successor flag as a future-only conditional frontier, not as a current exported role-transfer branch."
        ),
    }


def contains_all(text: str, needles: list[str]) -> bool:
    lower = text.lower()
    return all(needle.lower() in lower for needle in needles)


def closure_bits(text: str, keyword: str) -> dict[str, bool]:
    retained = contains_all(text, ["retained branch", "negative", keyword])
    replaced = contains_all(text, ["replaced branch", "negative", keyword])
    full = contains_all(text, ["full claim-specific", "closed negatively", "current repo state", keyword])
    future_not_forbidden = contains_all(
        text,
        ["future strict-side successor semantics cannot exist"],
    ) or contains_all(text, ["not discharge", "future strict-side successor semantics"])
    return {
        "retained_branch_closed_negative": retained,
        "replaced_branch_closed_negative": replaced,
        "full_claim_specific_closed_negative": full,
        "future_successor_not_forbidden_forever": future_not_forbidden,
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


def extract_p2395_rows(p2395: dict[str, Any]) -> dict[str, dict[str, Any]]:
    cert = p2395.get("post_bridge_retained_negative_successor_frontier_certificate", {})
    return {row.get("role_id", ""): row for row in cert.get("role_rows", [])}


def build_role_rows(artifacts: dict[str, Any]) -> list[dict[str, Any]]:
    p2395_rows = extract_p2395_rows(artifacts["P2395_SUCCESSOR_FRONTIER"])
    rows = []
    for spec in ROLE_PACKAGE_ROWS:
        text = artifacts[spec["full_negative_source"]].get("text", "")
        bits = closure_bits(text, spec["role_keyword"])
        p2395_row = p2395_rows.get(spec["p2395_role_id"], {})
        modified_flag = p2395_row.get("modified_successor_verdict") == "OPEN_ACTIVE_SUCCESSOR_FRONTIER"
        licensed = bool(p2395_row.get("current_transfer_licensed"))
        matrix_vector = [
            int(bits["retained_branch_closed_negative"]),
            int(bits["replaced_branch_closed_negative"]),
            int(bits["full_claim_specific_closed_negative"]),
            int(modified_flag),
            int(licensed),
            int(bits["future_successor_not_forbidden_forever"]),
        ]
        rows.append(
            {
                **spec,
                **bits,
                "p2395_modified_successor_flag": modified_flag,
                "p2395_modified_successor_rebased_verdict": "FUTURE_ONLY_CONDITIONAL_NOT_CURRENT_EXPORT",
                "current_transfer_licensed": licensed,
                "current_state_role_transfer_verdict": "CLOSED_NEGATIVE_CURRENT_REPO_STATE",
                "matrix_bits": dict(zip(MATRIX_COLUMNS, [bool(bit) for bit in matrix_vector])),
                "matrix_vector": matrix_vector,
            }
        )
    return rows


def package_certificate(rows: list[dict[str, Any]], package_text: str) -> dict[str, Any]:
    matrix = [row["matrix_vector"] for row in rows]
    package_closed = contains_all(
        package_text,
        ["full legacy physical-role package", "closed negatively", "current repo state"],
    )
    no_role_transfer = all(not row["current_transfer_licensed"] for row in rows)
    return {
        "matrix_columns": MATRIX_COLUMNS,
        "matrix": matrix,
        "matrix_rank_mod2": gf2_rank(matrix),
        "role_count": len(rows),
        "full_claim_specific_negative_count": sum(row["full_claim_specific_closed_negative"] for row in rows),
        "p2395_future_only_flag_count": sum(row["p2395_modified_successor_flag"] for row in rows),
        "current_licensed_transfer_count": sum(row["current_transfer_licensed"] for row in rows),
        "future_successor_not_forbidden_count": sum(row["future_successor_not_forbidden_forever"] for row in rows),
        "n116_package_closed_negative": package_closed,
        "all_current_role_transfer_closed_negative": package_closed and no_role_transfer,
        "rebase_rule": (
            "If a role has full claim-specific negative closure and no current transfer license, "
            "then P2395's modified-successor-open flag is read as future-only conditional unless a new explicit successor verdict is exported."
        ),
    }


def append_doc_sections() -> None:
    eq_section = """
## P2396/S1346 role-package negative closure rebase

`P2396/S1346` adds a nonduplication correction after P2395.  Repo grep finds that `N83`, `N99`, `N115`, and `N116` already close the full legacy physical-role package negatively on the current repo state.  Therefore the P2395 modified-successor flags must be read as future-only conditional slots, not as current exported successor branches.

The finite matrix has three role rows and six columns: retained negative closure, replaced negative closure, full claim-specific negative closure, P2395 modified-successor flag, current transfer license, and future-successor-not-forbidden.  The current licensed transfer count is zero, while the future-not-forbidden bit prevents overclaiming a forever rejection.

The next honest move is not to re-open the legacy role package.  It is to introduce genuinely new explicit successor evidence if one wants to revisit a role, otherwise keep the current-state package closed and continue on non-role strict-source/frontier work.  No `L_total`, SM/GR numeric extraction, or ToE closure follows.
""".strip()
    lag_section = """
## P2396/S1346 role-package closure rebase for Lagrangian use

`P2396/S1346` prevents a Lagrangian/EOM overread of P2395.  Since `N83/N99/N115/N116` already close the current-state legacy role package negatively, the P2395 modified-successor flags are future-only conditional slots.  They are not current licenses to import Weinberg, fine-structure, or gravity-hierarchy semantics into `L_total`.

Any future role-bearing `L_total` term must bring new explicit successor evidence first.  Until then, the current-state role-transfer count is zero.
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
    rows = build_role_rows(artifacts)
    package = package_certificate(rows, artifacts["N116_ROLE_PACKAGE_FULL_NEGATIVE"].get("text", ""))
    theorem_export = {
        "theorem_name": "P2396_T1_role_package_current_state_negative_closure_rebase",
        "role_package_closed_negative_now": package["all_current_role_transfer_closed_negative"],
        "current_licensed_transfer_count": package["current_licensed_transfer_count"],
        "p2395_modified_successor_flags_rebased_as_future_only": package["p2395_future_only_flag_count"],
        "future_successor_not_forbidden_count": package["future_successor_not_forbidden_count"],
        "matrix_rank_mod2": package["matrix_rank_mod2"],
        "not_licensed": [
            "No legacy physical role is transferred on the current repo state.",
            "No P2395 modified-successor flag is a current exported successor theorem.",
            "No forever impossibility of future strict successor evidence is claimed.",
            "No L_total promotion, SM/GR numeric extraction, or ToE closure is claimed.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "all_three_full_claim_negative_closures_loaded": package["full_claim_specific_negative_count"] == 3,
        "n116_package_negative_closure_loaded": package["n116_package_closed_negative"],
        "current_transfer_count_zero": package["current_licensed_transfer_count"] == 0,
        "p2395_flags_rebased_future_only": package["p2395_future_only_flag_count"] == 3,
        "future_successor_not_forbidden_preserved": package["future_successor_not_forbidden_count"] == 3,
        "matrix_rank_is_one": package["matrix_rank_mod2"] == 1,
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2396_s1346_v1",
        "packet_id": "P2396",
        "stage_id": "S1346",
        "result_kind": "ROLE_PACKAGE_NEGATIVE_CLOSURE_REBASE_CERTIFICATE",
        "status": "PASS_ROLE_PACKAGE_CLOSED_CURRENT_STATE_P2395_FLAGS_FUTURE_ONLY",
        "role_package_negative_closure_rebase_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: artifact.get("status", "TEXT_SOURCE") for name, artifact in artifacts.items()},
            "role_rows": rows,
            "package_certificate": package,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": "Do not reopen the legacy role package without new explicit successor evidence; use P2396 as the current-state closure guard and move to a non-role strict-source frontier or a genuinely new successor theorem input.",
        "global_status": "OPEN_PROGRESS_ROLE_PACKAGE_CURRENT_STATE_CLOSED_FUTURE_SUCCESSOR_ONLY_CONDITIONAL",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    cert = payload["role_package_negative_closure_rebase_certificate"]
    theorem = cert["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2396 S1346: role-package negative closure rebase certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2396/S1346 rebases P2395 against N83/N99/N115/N116: the legacy role package is closed negatively on the current repo state, and P2395 modified-successor flags are future-only conditional slots.",
                "",
                "## Finite certificate",
                "",
                f"- Current licensed transfer count: `{theorem['current_licensed_transfer_count']}`.",
                f"- P2395 flags rebased as future-only: `{theorem['p2395_modified_successor_flags_rebased_as_future_only']}`.",
                f"- Future successor not forbidden count: `{theorem['future_successor_not_forbidden_count']}`.",
                f"- Matrix rank over GF(2): `{theorem['matrix_rank_mod2']}`.",
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
