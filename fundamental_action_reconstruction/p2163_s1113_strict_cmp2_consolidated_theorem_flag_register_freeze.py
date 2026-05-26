#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2162 = GEN / "p2162_s1112_strict_cmp2_end_to_end_refresh_and_theorem_flag_ledger_snapshot.json"
OUT = GEN / "p2163_s1113_strict_cmp2_consolidated_theorem_flag_register_freeze.json"
MD = GEN / "p2163_s1113_strict_cmp2_consolidated_theorem_flag_register_freeze.md"


def load(p: Path) -> dict[str, Any]:
    if not p.exists():
        return {"_missing": str(p), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2162 = load(IN_2162)
    block = p2162.get("end_to_end_refresh_and_theorem_flag_ledger_snapshot", {}) or {}
    source_truth = block.get("source_of_truth_flags", {}) or {}
    rows = block.get("ledger_rows", []) or []

    consolidated = {
        r.get("packet"): {
            "result_kind": r.get("result_kind"),
            "full_d3_covariance_transport_proven": r.get("full_d3_covariance_transport_proven"),
            "c3_theorem_proven": r.get("c3_theorem_proven"),
            "status": r.get("status"),
        }
        for r in rows
    }

    snapshot_consistent = bool(block.get("snapshot_consistent", False))
    result_kind = (
        "PASS_STRICT_CMP2_CONSOLIDATED_THEOREM_FLAG_REGISTER_FREEZE"
        if snapshot_consistent
        else "OPEN_STRICT_CMP2_CONSOLIDATED_THEOREM_FLAG_REGISTER_FREEZE_WITH_DRIFT"
    )

    payload = {
        "schema_version": "p2163_s1113_v1",
        "packet_id": "P2163",
        "stage_id": "S1113",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": result_kind,
        "consolidated_theorem_flag_register_freeze": {
            "source_snapshot": str(IN_2162.relative_to(ROOT)),
            "source_of_truth_flags": source_truth,
            "packet_register": consolidated,
            "n_packets": len(consolidated),
            "snapshot_consistent": snapshot_consistent,
            "freeze_label": "STRICT_CMP2_D3_C3_FLAG_BASELINE_V1",
            "scope_limit": "consolidated register only; not global ToE closure",
        },
        "recommended_next_honest_step": {
            "id": "P2164_candidate",
            "goal": "start next theorem blocker lane (e.g., Cutkosky/unitarity or selector/QW-2191) from this frozen baseline",
        },
        "gatekeeper_checks": {
            "register_freeze_exported": True,
            "snapshot_consistent": snapshot_consistent,
            "no_toe_closure_claimed": True,
            "full_d3_covariance_transport_proven": bool(source_truth.get("full_d3_covariance_transport_proven", False)),
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": bool(source_truth.get("c3_theorem_proven", False)),
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2163 S1113: strict CMP2 consolidated theorem-flag register freeze",
                "",
                f"- Result kind: `{result_kind}`",
                f"- n_packets: `{len(consolidated)}`",
                f"- snapshot_consistent: `{snapshot_consistent}`",
                "",
                "No global ToE closure claim is made.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
