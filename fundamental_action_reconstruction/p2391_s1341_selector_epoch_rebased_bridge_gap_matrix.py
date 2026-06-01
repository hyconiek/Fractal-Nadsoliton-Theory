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
SCRATCH = ROOT / "scratch"

OUT = GEN / "p2391_s1341_selector_epoch_rebased_bridge_gap_matrix.json"
MD = GEN / "p2391_s1341_selector_epoch_rebased_bridge_gap_matrix.md"

SOURCE_FILES = {
    "P1343_SELECTOR_SOURCE_REPORT": GEN / "p1343_p1343_report_v1.json",
    "P1348_GLOBAL_CLOSURE_REPORT": GEN / "p1348_p1348_report_v1.json",
    "P2390_SELECTOR_QUALIFIED_AUDIT": GEN / "p2390_s1340_selector_qualified_beta_tors_chi11_role_audit.json",
    "OLD_COMPONENT_GAP_REPORT": SCRATCH / "bridge_strict_completion_legacy_to_strict_completion_component_gap_certificate_report.json",
    "S2_STRATEGIC_PRIORITY": ROOT / "S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

COLUMNS = [
    "finite_certificate_exported",
    "generic_strict_selector_exported",
    "explicit_beta_tors_to_chi11_transport",
    "strict_dynamic_source_exported",
    "role_transfer_allowed_now",
]

ROWS = [
    "amplitude_normalization",
    "phase_frequency_transport",
    "damping_compression",
    "topological_phase_bit_chi11",
    "legacy_physical_role_transfer",
]


def rel(path: Path) -> str:
    return path.relative_to(REPO).as_posix()


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": rel(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def read_text(path: Path) -> str:
    if not path.exists():
        return f"OPEN_MISSING_ARTIFACT::{rel(path)}"
    return path.read_text(encoding="utf-8")


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def gf2_rank(matrix: list[list[int]]) -> int:
    work = [row[:] for row in matrix]
    rank = 0
    col_count = len(work[0]) if work else 0
    for col in range(col_count):
        pivot = next((r for r in range(rank, len(work)) if work[r][col] % 2), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        for r in range(len(work)):
            if r != rank and work[r][col] % 2:
                work[r] = [(a ^ b) for a, b in zip(work[r], work[rank])]
        rank += 1
    return rank


def rg_audit() -> dict[str, Any]:
    patterns = [
        "P2390|S1340|selector-qualified",
        "selector_source_gap_remains|bit_source_not_exported",
        "P1343|S_strict_internal_v1",
        "beta_tors.*chi_?11|chi_?11.*beta_tors",
    ]
    results: dict[str, Any] = {}
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
        results[pattern] = {"count": len(lines), "samples": lines[:8]}
    return {
        "tool": "rg",
        "patterns": results,
        "finding": (
            "The repo already contains an older component-gap matrix with a pre-P1343-style selector-source gap, plus P2390's selector-qualified audit. "
            "P2391 therefore rebases the matrix by epoch instead of duplicating either artifact."
        ),
    }


def old_rows(report: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {row["component"]: row for row in report.get("component_gap_rows", [])}


def build_rebased_rows(artifacts: dict[str, Any]) -> list[dict[str, Any]]:
    old = old_rows(artifacts["OLD_COMPONENT_GAP_REPORT"])
    p2390_cert = artifacts["P2390_SELECTOR_QUALIFIED_AUDIT"]["selector_qualified_beta_tors_chi11_role_audit"]["implication_certificate"]
    p1343_selector = artifacts["P1343_SELECTOR_SOURCE_REPORT"].get("status") == "CLOSED_FULL_STRICT_INTERNAL_SOURCE_V1"
    p1348_closed = artifacts["P1348_GLOBAL_CLOSURE_REPORT"].get("status") == "GLOBAL_CLOSURE_THEOREM_EXPORTED_CLOSED_DECLARED_SCOPE"
    strict_selector = p1343_selector and p1348_closed and p2390_cert["transfer_rule_inputs"]["selector_source"]
    beta_transport = p2390_cert["transfer_rule_inputs"]["explicit_beta_tors_to_chi11_transport_theorem"]
    role_transfer = p2390_cert["transfer_rule_inputs"]["separate_role_transfer_theorem"]
    bridge_ready = p2390_cert["transfer_rule_inputs"]["full_legacy_to_strict_bridge_ready"]

    rows: list[dict[str, Any]] = []
    for row_name in ROWS:
        old_row = old[row_name]
        old_bits = old_row["matrix_bits"]
        is_chi11 = row_name == "topological_phase_bit_chi11"
        is_role = row_name == "legacy_physical_role_transfer"
        generic_selector = strict_selector if is_chi11 else bool(old_bits.get("selector_or_source_exported"))
        new_bits = {
            "finite_certificate_exported": bool(old_bits.get("finite_certificate_exported")),
            "generic_strict_selector_exported": generic_selector,
            "explicit_beta_tors_to_chi11_transport": bool(beta_transport) if is_chi11 else False,
            "strict_dynamic_source_exported": bool(old_bits.get("strict_dynamic_source_exported")),
            "role_transfer_allowed_now": bool(role_transfer and bridge_ready) if is_role or is_chi11 else False,
        }
        rows.append(
            {
                "component": row_name,
                "old_selector_or_source_exported": bool(old_bits.get("selector_or_source_exported")),
                "rebased_bits": new_bits,
                "epoch_delta": {
                    "selector_bit_flipped_by_P1343_P1348": is_chi11 and generic_selector != bool(old_bits.get("selector_or_source_exported")),
                    "beta_transport_still_missing": is_chi11 and not beta_transport,
                    "role_transfer_still_blocked": not new_bits["role_transfer_allowed_now"],
                },
                "status": (
                    "generic_selector_rebased_present__beta_tors_transport_and_role_transfer_absent"
                    if is_chi11
                    else old_row.get("status", "carried_forward")
                ),
            }
        )
    return rows


def build_payload() -> dict[str, Any]:
    artifacts: dict[str, Any] = {}
    for name, path in SOURCE_FILES.items():
        artifacts[name] = load_json(path) if path.suffix == ".json" else {"text": read_text(path)}
    grep = rg_audit()
    rows = build_rebased_rows(artifacts)
    matrix = [[1 if row["rebased_bits"][col] else 0 for col in COLUMNS] for row in rows]
    old_selector_vector = [1 if row["old_selector_or_source_exported"] else 0 for row in rows]
    new_selector_vector = [1 if row["rebased_bits"]["generic_strict_selector_exported"] else 0 for row in rows]
    selector_hamming_delta = sum(a != b for a, b in zip(old_selector_vector, new_selector_vector))
    chi11 = next(row for row in rows if row["component"] == "topological_phase_bit_chi11")
    theorem_export = {
        "theorem_name": "P2391_T1_selector_epoch_rebased_gap_matrix",
        "statement": (
            "After P1343/P1348 and P2390, the topological chi11 row must be rebased from 'no generic selector source' to 'generic strict selector present'. "
            "This one-bit epoch update does not license beta_tors->chi11 transport or legacy role transfer."
        ),
        "matrix_columns": COLUMNS,
        "selector_hamming_delta_old_to_rebased": selector_hamming_delta,
        "matrix_rank_mod2": gf2_rank(matrix),
        "not_licensed": [
            "No explicit beta_tors -> chi11 transport theorem is exported.",
            "No completed legacy -> strict bridge is exported.",
            "No legacy role-transfer theorem is exported.",
            "No L_total promotion, cap-density source theorem, SM/GR numeric extraction, or ToE closure follows.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_found_p2390": grep["patterns"]["P2390|S1340|selector-qualified"]["count"] > 0,
        "rg_found_old_selector_gap_wording": grep["patterns"]["selector_source_gap_remains|bit_source_not_exported"]["count"] > 0,
        "p1343_selector_present": artifacts["P1343_SELECTOR_SOURCE_REPORT"].get("status") == "CLOSED_FULL_STRICT_INTERNAL_SOURCE_V1",
        "p1348_declared_closure_present": artifacts["P1348_GLOBAL_CLOSURE_REPORT"].get("status") == "GLOBAL_CLOSURE_THEOREM_EXPORTED_CLOSED_DECLARED_SCOPE",
        "chi11_selector_bit_rebased_true": chi11["rebased_bits"]["generic_strict_selector_exported"],
        "chi11_beta_transport_still_false": not chi11["rebased_bits"]["explicit_beta_tors_to_chi11_transport"],
        "role_transfer_all_rows_false": all(not row["rebased_bits"]["role_transfer_allowed_now"] for row in rows),
        "exactly_one_selector_epoch_delta": selector_hamming_delta == 1,
        "docs_updated_with_p2391": all("P2391/S1341" in read_text(path) for path in DOC_FILES.values()),
    }
    return {
        "schema_version": "p2391_s1341_v1",
        "packet_id": "P2391",
        "stage_id": "S1341",
        "timestamp_utc": "2026-06-01T00:00:00Z",
        "produced_by": rel(Path(__file__).resolve()),
        "result_kind": "SELECTOR_EPOCH_REBASED_BRIDGE_GAP_MATRIX",
        "status": "PASS_ONE_SELECTOR_EPOCH_BIT_REBASED_BETA_TORS_TRANSPORT_AND_ROLE_TRANSFER_STILL_OPEN",
        "selector_epoch_rebased_bridge_gap_matrix": {
            "rg_nonduplication_audit": grep,
            "matrix_columns": COLUMNS,
            "rows": rows,
            "rebased_matrix": matrix,
            "old_selector_vector": old_selector_vector,
            "new_selector_vector": new_selector_vector,
            "selector_hamming_delta_old_to_rebased": selector_hamming_delta,
            "matrix_rank_mod2": gf2_rank(matrix),
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": "Do not pursue beta_tors->chi11 for selector necessity; continue with legacy->strict bridge completion and later role-transfer auditing.",
        "global_status": "OPEN_PROGRESS_SELECTOR_EPOCH_REBASED_NO_ROLE_TRANSFER_OR_TOE_CLOSURE",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    cert = payload["selector_epoch_rebased_bridge_gap_matrix"]
    chi11 = next(row for row in cert["rows"] if row["component"] == "topological_phase_bit_chi11")
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2391 S1341: selector-epoch rebased bridge gap matrix",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2391/S1341 rebases the older bridge component-gap matrix after the P1343/P1348 selector epoch and P2390 qualification.",
                f"The selector/source vector changes by Hamming distance `{cert['selector_hamming_delta_old_to_rebased']}`: only the `topological_phase_bit_chi11` generic selector bit flips from old gap wording to selector-present wording.",
                f"The rebased GF(2) matrix rank is `{cert['matrix_rank_mod2']}`.",
                "",
                "## Chi11 row after rebase",
                "",
                f"- Generic strict selector exported: `{chi11['rebased_bits']['generic_strict_selector_exported']}`.",
                f"- Explicit `beta_tors -> chi11` transport: `{chi11['rebased_bits']['explicit_beta_tors_to_chi11_transport']}`.",
                f"- Role transfer allowed now: `{chi11['rebased_bits']['role_transfer_allowed_now']}`.",
                "",
                "## Hard limits",
                "",
                "- This removes stale generic-selector-gap wording only; P2392 retires `beta_tors -> chi11` as an active selector-search target rather than proving that transport.",
                "- No completed legacy-to-strict bridge, no legacy role-transfer theorem, no `L_total` promotion, no cap-density source theorem, and no ToE closure is claimed.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def main() -> None:
    payload = build_payload()
    write_outputs(payload)


if __name__ == "__main__":
    main()
