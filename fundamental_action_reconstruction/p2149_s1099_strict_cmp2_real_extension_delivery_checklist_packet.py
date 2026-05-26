#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2148 = GEN / "p2148_s1098_strict_cmp2_external_data_delivery_blocker_report.json"
OUT = GEN / "p2149_s1099_strict_cmp2_real_extension_delivery_checklist_packet.json"
MD = GEN / "p2149_s1099_strict_cmp2_real_extension_delivery_checklist_packet.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path)}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2148 = load(IN_2148)
    blk = (p2148.get("external_data_delivery_blocker_report", {}) or {})
    blocked = blk.get("blocked_stages", []) or []

    checklist = [
        "prepare generated/cmp2_backend_rows_extension_v1.json",
        "ensure each row has keys: cmp_bin_index, backend_s, posterior_weights_backend_evidence, posterior_predictive_covered",
        "provide at least 8 valid new rows",
        "rerun p2147 checkpoint",
        "confirm PASS on p2132/p2133/p2134/p2135 before interpretation",
    ]

    payload = {
        "schema_version": "p2149_s1099_v1",
        "packet_id": "P2149",
        "stage_id": "S1099",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_CMP2_REAL_EXTENSION_DELIVERY_CHECKLIST_PACKET_WITH_TRACE",
        "real_extension_delivery_checklist": {
            "source_blocker_report": str(IN_2148.relative_to(ROOT)),
            "blocked_stages": blocked,
            "checklist": checklist,
            "scope_limit": "delivery checklist only; not theorem-grade closure",
        },
        "recommended_next_honest_step": {
            "id": "P2150_candidate",
            "goal": "after checklist completion rerun P2147 and then P2141/P2146 for refreshed real-data flags",
        },
        "gatekeeper_checks": {
            "checklist_exported": True,
            "no_toe_closure_claimed": True,
            "full_d3_covariance_transport_proven": False,
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": False,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2149 S1099: strict CMP2 real extension delivery checklist packet\n\n"
        f"- Result kind: `{payload['result_kind']}`\n"
        f"- Blocked stages: `{blocked}`\n\n"
        "No theorem-grade global closure claim is made.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
