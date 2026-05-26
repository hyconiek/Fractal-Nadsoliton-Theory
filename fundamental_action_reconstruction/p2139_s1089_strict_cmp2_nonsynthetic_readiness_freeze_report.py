#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2138 = GEN / "p2138_s1088_strict_cmp2_nonsynthetic_rerun_orchestrator.json"
OUT = GEN / "p2139_s1089_strict_cmp2_nonsynthetic_readiness_freeze_report.json"
MD = GEN / "p2139_s1089_strict_cmp2_nonsynthetic_readiness_freeze_report.md"

SCHEMA_VERSION = "p2139_s1089_v1"
TIMESTAMP_UTC = "2026-05-26T00:00:00+00:00"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path)}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2138 = load(IN_2138)

    orch = p2138.get("nonsynthetic_rerun_orchestrator", {}) or {}
    rk = (orch.get("post_run_result_kinds", {}) or {})

    pass_2132 = rk.get("p2132") == "PASS_STRICT_CMP2_REAL_SUPPORT_ACQUISITION_GATE_WITH_TRACE"
    pass_2133 = rk.get("p2133") == "PASS_STRICT_CMP2_REAL_EXTENSION_MERGE_CONTRACT_WITH_TRACE"
    pass_2134 = rk.get("p2134") == "PASS_STRICT_CMP2_NONSYNTHETIC_RERUN_COMPARISON_AUDIT_WITH_TRACE"
    pass_2135 = rk.get("p2135") == "PASS_STRICT_CMP2_MERGED_REAL_BLOCK_VARIANT_STABILITY_AUDIT_WITH_TRACE"
    fully_unblocked = pass_2132 and pass_2133 and pass_2134 and pass_2135

    missing = []
    if not pass_2132:
        missing.append("P2132")
    if not pass_2133:
        missing.append("P2133")
    if not pass_2134:
        missing.append("P2134")
    if not pass_2135:
        missing.append("P2135")

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2139",
        "stage_id": "S1089",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_CMP2_NONSYNTHETIC_READINESS_FREEZE_REPORT_WITH_TRACE" if fully_unblocked else "OPEN_STRICT_CMP2_NONSYNTHETIC_READINESS_FREEZE_REPORT_BLOCKED",
        "nonsynthetic_readiness_freeze_report": {
            "source_orchestrator": str(IN_2138.relative_to(ROOT)),
            "post_run_result_kinds": rk,
            "fully_unblocked": fully_unblocked,
            "missing_or_blocked_stages": missing,
            "scope_limit": "readiness freeze/readout only; not theorem-grade closure",
        },
        "recommended_next_honest_step": {
            "id": "P2140_candidate",
            "goal": "if blocked: deliver real extension rows and rerun P2138; if unblocked: lock non-synthetic CI inflation report and move to next formal obstruction family",
        },
        "gatekeeper_checks": {
            "all_nonsynthetic_gates_passed": fully_unblocked,
            "full_d3_covariance_transport_proven": False,
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": False,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2139 S1089: strict CMP2 non-synthetic readiness freeze report",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Fully unblocked: `{fully_unblocked}`",
            "",
            "This stage freezes readiness status after P2138 orchestration and lists remaining blockers.",
            "No theorem-grade global closure claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
