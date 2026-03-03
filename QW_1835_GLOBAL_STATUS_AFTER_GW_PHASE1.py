#!/usr/bin/env python3
"""
QW-1835: Global empirical status update after GW redesign phase-1.

Compares pre-update global status with new GW phase-1 gate.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1835_global_status_after_gw_phase1.json"
OUT_MD = ROOT / "RAPORT_QW1835_GLOBAL_STATUS_AFTER_GW_PHASE1.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    d1824 = load("report_qw1824_quantile_gated_readiness.json")
    d1825 = load("report_qw1825_quantile_protocol_transfer_gate.json")
    d1830 = load("report_qw1830_global_empirical_status_gate.json")
    d1834 = load("report_qw1834_gw_redesign_phase1_gate.json")

    checks: List[Dict[str, object]] = []

    pta_ready = d1824.get("hard_gate") == "PASS"
    checks.append(
        {
            "domain": "PTA quantile branch",
            "score": 1.0 if pta_ready else 0.0,
            "status": "PASS" if pta_ready else "FAIL",
            "note": d1824.get("readiness", ""),
        }
    )

    transfer = d1825.get("hard_gate", "FAIL")
    score_transfer = 1.0 if transfer == "PASS" else (0.6 if transfer == "PARTIAL" else 0.0)
    checks.append(
        {
            "domain": "Cross-domain transfer",
            "score": score_transfer,
            "status": transfer,
            "note": d1825.get("readiness", ""),
        }
    )

    gw_hard = d1834.get("hard_gate", "FAIL")
    score_gw = 1.0 if gw_hard == "PASS" else (0.6 if gw_hard == "PARTIAL" else 0.0)
    checks.append(
        {
            "domain": "GW redesign phase-1",
            "score": score_gw,
            "status": gw_hard,
            "note": d1834.get("readiness", ""),
        }
    )

    new_score = float(sum(float(c["score"]) for c in checks) / len(checks))
    old_score = float(d1830.get("global_score", 0.0))
    delta = new_score - old_score

    if pta_ready and transfer in ("PASS", "PARTIAL") and gw_hard == "PASS":
        readiness = "GLOBAL_READY_FOR_JOINT_CONFIRMATION"
        hard_gate = "PASS"
        recommendation = "RUN_JOINT_PTA_GW_CONFIRMATORY_PROTOCOL"
    elif pta_ready and transfer in ("PASS", "PARTIAL") and gw_hard == "PARTIAL":
        readiness = "GLOBAL_PARTIAL_WITH_GW_PHASE2_REQUIRED"
        hard_gate = "PARTIAL"
        recommendation = "LAUNCH_GW_PHASE2_STRUCTURAL_REPARAM_THEN_RETEST_GLOBAL_GATE"
    else:
        readiness = "GLOBAL_NOT_READY"
        hard_gate = "FAIL"
        recommendation = "CONTINUE_METHOD_REPAIR"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "global_score_previous_qw1830": old_score,
        "global_score_updated": new_score,
        "delta_global_score": delta,
        "hard_gate": hard_gate,
        "readiness": readiness,
        "recommendation": recommendation,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1835: GLOBAL STATUS AFTER GW PHASE-1",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Global score (QW-1830 -> now): {old_score:.3f} -> {new_score:.3f}",
        f"- Delta score: {delta:.3f}",
        f"- Hard gate: **{hard_gate}**",
        f"- Readiness: **{readiness}**",
        f"- Recommendation: **{recommendation}**",
        "",
        "## Checks",
    ]
    for c in checks:
        lines.append(f"- {c['domain']}: {c['status']} | score={c['score']:.3f} | note={c['note']}")

    lines.extend(["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1835] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1835] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
