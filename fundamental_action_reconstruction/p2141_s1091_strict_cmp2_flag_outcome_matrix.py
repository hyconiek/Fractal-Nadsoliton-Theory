#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2141_s1091_strict_cmp2_flag_outcome_matrix.json"
MD = GEN / "p2141_s1091_strict_cmp2_flag_outcome_matrix.md"

SCHEMA_VERSION = "p2141_s1091_v1"
TIMESTAMP_UTC = "2026-05-26T00:00:00+00:00"

STAGES = [
    "p2132_s1082_strict_cmp2_real_support_acquisition_gate",
    "p2133_s1083_strict_cmp2_real_extension_merge_contract",
    "p2134_s1084_strict_cmp2_nonsynthetic_rerun_comparison_audit",
    "p2135_s1085_strict_cmp2_merged_real_block_variant_stability_audit",
    "p2138_s1088_strict_cmp2_nonsynthetic_rerun_orchestrator",
    "p2139_s1089_strict_cmp2_nonsynthetic_readiness_freeze_report",
    "p2140_s1090_strict_cmp2_blocked_path_resolution_router",
]


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)

    matrix = {}
    for stem in STAGES:
        payload = load(GEN / f"{stem}.json")
        matrix[stem] = {
            "result_kind": payload.get("result_kind", "OPEN_UNKNOWN"),
            "status": payload.get("status", "OPEN_UNKNOWN"),
            "gatekeeper_checks": payload.get("gatekeeper_checks", {}),
        }

    pass_count = sum(1 for v in matrix.values() if str(v["result_kind"]).startswith("PASS_"))
    open_count = sum(1 for v in matrix.values() if str(v["result_kind"]).startswith("OPEN_"))

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2141",
        "stage_id": "S1091",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_CMP2_FLAG_OUTCOME_MATRIX_WITH_TRACE",
        "flag_outcome_matrix": {
            "tracked_stages": STAGES,
            "matrix": matrix,
            "pass_count": pass_count,
            "open_count": open_count,
            "scope_limit": "flag/readout matrix only; not theorem-grade closure",
        },
        "recommended_next_honest_step": {
            "id": "P2142_candidate",
            "goal": "resolve OPEN_* flags in P2132/P2133 first by delivering real extension rows, then rerun P2138-P2141",
        },
        "gatekeeper_checks": {
            "matrix_exported": True,
            "no_toe_closure_claimed": True,
            "full_d3_covariance_transport_proven": False,
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": False,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    lines = [
        "# P2141 S1091: strict CMP2 flag outcome matrix",
        "",
        f"- Status: `{payload['status']}`",
        f"- Result kind: `{payload['result_kind']}`",
        f"- PASS count: `{pass_count}`",
        f"- OPEN count: `{open_count}`",
        "",
        "## Key flag readout",
    ]
    for stem, v in matrix.items():
        lines.append(f"- `{stem}` -> `{v['result_kind']}`")
    lines.extend([
        "",
        "No theorem-grade global closure claim is made.",
    ])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
