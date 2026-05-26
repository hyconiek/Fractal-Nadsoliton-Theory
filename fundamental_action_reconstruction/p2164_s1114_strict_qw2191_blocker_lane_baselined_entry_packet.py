#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2163 = GEN / "p2163_s1113_strict_cmp2_consolidated_theorem_flag_register_freeze.json"
OUT = GEN / "p2164_s1114_strict_qw2191_blocker_lane_baselined_entry_packet.json"
MD = GEN / "p2164_s1114_strict_qw2191_blocker_lane_baselined_entry_packet.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2163 = load(IN_2163)
    b = p2163.get("consolidated_theorem_flag_register_freeze", {}) or {}

    freeze_label = b.get("freeze_label")
    snapshot_consistent = bool(b.get("snapshot_consistent", False))
    source_flags = b.get("source_of_truth_flags", {}) or {}

    baseline_ready = snapshot_consistent and freeze_label == "STRICT_CMP2_D3_C3_FLAG_BASELINE_V1"

    result_kind = (
        "PASS_STRICT_QW2191_BLOCKER_LANE_BASELINED_ENTRY_PACKET"
        if baseline_ready
        else "OPEN_STRICT_QW2191_BLOCKER_LANE_BASELINED_ENTRY_PACKET_BLOCKED"
    )

    payload = {
        "schema_version": "p2164_s1114_v1",
        "packet_id": "P2164",
        "stage_id": "S1114",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": result_kind,
        "strict_qw2191_blocker_lane_baselined_entry_packet": {
            "source_baseline": str(IN_2163.relative_to(ROOT)),
            "baseline_freeze_label": freeze_label,
            "baseline_snapshot_consistent": snapshot_consistent,
            "d3_c3_source_flags": source_flags,
            "selected_next_blocker_lane": "QW-2191_selector_obstruction",
            "entry_contract": {
                "respect_strict_lane_only": True,
                "no_selector_closure_claim": True,
                "no_l5_l12_cyclic_expansion": True,
                "baseline_required": baseline_ready,
            },
            "scope_limit": "baselined entry packet for next blocker lane only; no global ToE closure",
        },
        "recommended_next_honest_step": {
            "id": "P2165_candidate",
            "goal": "export strict selector-premise hypothesis lattice and noncyclic admissibility checker for QW-2191 lane",
        },
        "gatekeeper_checks": {
            "baselined_entry_exported": True,
            "baseline_snapshot_consistent": snapshot_consistent,
            "no_toe_closure_claimed": True,
            "full_d3_covariance_transport_proven": bool(source_flags.get("full_d3_covariance_transport_proven", False)),
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": bool(source_flags.get("c3_theorem_proven", False)),
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2164 S1114: strict QW-2191 blocker-lane baselined entry packet",
                "",
                f"- Result kind: `{result_kind}`",
                f"- Baseline ready: `{baseline_ready}`",
                "- Selected blocker lane: `QW-2191_selector_obstruction`",
                "",
                "No global ToE closure claim is made.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
