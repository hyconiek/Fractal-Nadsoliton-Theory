#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2201 = GEN / "p2201_s1151_strict_cutkosky_shared_majorant_frw_bianchi_certificate.json"
IN_2025 = GEN / "p2025_s975_strict_cutkosky_same_scheme_cohomology_amplitude_bridge_seed.json"
OUT = GEN / "p2202_s1152_strict_cutkosky_frw_bianchi_bound_margin_vs_task3_obstruction_ledger.json"
MD = GEN / "p2202_s1152_strict_cutkosky_frw_bianchi_bound_margin_vs_task3_obstruction_ledger.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2201 = load(IN_2201)
    p2025 = load(IN_2025)

    cert = p2201.get("strict_cutkosky_shared_majorant_frw_bianchi_certificate", {}) or {}
    summary = cert.get("margin_summary", {}) or {}
    rows = cert.get("rows", []) or []

    task_rows = p2025.get("toe_closure_gaps_7tasks", []) or []
    task3 = next((r for r in task_rows if int(r.get("task_id", -1)) == 3), {})

    min_margin = float(summary.get("min_margin", 0.0) or 0.0)
    spread = float(summary.get("margin_spread", 0.0) or 0.0)
    all_bounds = bool(summary.get("all_bounds_respected", False))

    evidence_grade = (
        "quantitative_support_for_nonnegative_bound_margin"
        if all_bounds and min_margin >= 0.0
        else "insufficient_bound_margin_support"
    )

    obstruction_compatible = str(task3.get("status", "OPEN_OBSTRUCTION_WITH_TRACE")) == "OPEN_OBSTRUCTION_WITH_TRACE"

    ledger = {
        "ledger_id": "STRICT_FRW_BIANCHI_BOUND_MARGIN_VS_TASK3_OBSTRUCTION_LEDGER_V1",
        "source_packets": [str(IN_2201.relative_to(ROOT)), str(IN_2025.relative_to(ROOT))],
        "task3_row": {
            "task_id": task3.get("task_id", 3),
            "name": task3.get("name", "Background-independence global transport"),
            "status": task3.get("status", "OPEN_OBSTRUCTION_WITH_TRACE"),
            "missing": task3.get("missing", []),
        },
        "bound_margin_evidence": {
            "rows": rows,
            "min_margin": min_margin,
            "margin_spread": spread,
            "all_bounds_respected": all_bounds,
            "evidence_grade": evidence_grade,
        },
        "cross_verdict": {
            "obstruction_status_retained": obstruction_compatible,
            "does_not_discharge_task3": True,
            "reason": "Bound margins strengthen quantitative control but do not provide theorem-grade FRW<->Bianchi transport closure.",
        },
    }

    payload = {
        "schema_version": "p2202_s1152_v1",
        "packet_id": "P2202",
        "stage_id": "S1152",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_FRW_BIANCHI_BOUND_MARGIN_VS_TASK3_OBSTRUCTION_LEDGER",
        "strict_cutkosky_frw_bianchi_bound_margin_vs_task3_obstruction_ledger": ledger,
        "recommended_next_honest_step": {
            "id": "P2203_candidate",
            "goal": "add executable FRW<->Bianchi transport residual map using shared majorant controls and report theorem-gap deltas",
        },
        "gatekeeper_checks": {
            "obstruction_ledger_exported": True,
            "task3_status_preserved_open": obstruction_compatible,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_cutkosky_closure_proven": False,
            "full_d3_covariance_transport_proven": False,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2202 S1152: FRW/Bianchi bound-margin vs Task-3 obstruction ledger",
            "",
            f"- Task-3 status: `{task3.get('status', 'OPEN_OBSTRUCTION_WITH_TRACE')}`",
            f"- Min margin: `{min_margin:.12e}`",
            f"- Margin spread: `{spread:.12e}`",
            f"- All bounds respected: `{all_bounds}`",
            "",
            "This is an obstruction-ledger crosswalk only, not theorem-grade transport closure.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
