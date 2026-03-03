#!/usr/bin/env python3
"""
QW-1830: Global empirical status gate after quantile transition.

Aggregates:
- QW-1824 (PTA sequence readiness under quantile gate),
- QW-1825 (cross-domain transfer gate),
- QW-1827 (GW redesign gate),
- QW-1829 (GW near-target feasibility).
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1830_global_empirical_status_gate.json"
OUT_MD = ROOT / "RAPORT_QW1830_GLOBAL_EMPIRICAL_STATUS_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    d1824 = load("report_qw1824_quantile_gated_readiness.json")
    d1825 = load("report_qw1825_quantile_protocol_transfer_gate.json")
    d1827 = load("report_qw1827_gw_redesign_gate.json")
    d1829 = load("report_qw1829_gw_near_target_feasibility.json")

    checks: List[Dict[str, object]] = []

    pta_ready = d1824.get("hard_gate") == "PASS"
    checks.append(
        {
            "domain": "PTA quantile-gated branch",
            "score": 1.0 if pta_ready else 0.0,
            "status": "PASS" if pta_ready else "FAIL",
            "note": d1824.get("readiness", ""),
        }
    )

    transfer_partial = d1825.get("hard_gate") == "PARTIAL"
    transfer_pass = d1825.get("hard_gate") == "PASS"
    transfer_score = 1.0 if transfer_pass else (0.6 if transfer_partial else 0.0)
    checks.append(
        {
            "domain": "Cross-domain transfer status",
            "score": transfer_score,
            "status": "PASS" if transfer_pass else ("PARTIAL" if transfer_partial else "FAIL"),
            "note": d1825.get("readiness", ""),
        }
    )

    gw_redesign_required = d1827.get("readiness") == "GW_BRANCH_REDESIGN_REQUIRED"
    gw_score = float(d1827.get("global_score", 0.0))
    checks.append(
        {
            "domain": "GW branch readiness",
            "score": gw_score,
            "status": "FAIL" if gw_redesign_required else "PASS",
            "note": d1827.get("readiness", ""),
        }
    )

    structural_redesign = d1829.get("verdict") == "GW_NEAR_TARGET_REQUIRES_STRUCTURAL_REDESIGN"
    checks.append(
        {
            "domain": "GW near-target feasibility",
            "score": 0.0 if structural_redesign else 0.8,
            "status": "FAIL" if structural_redesign else "PASS",
            "note": d1829.get("verdict", ""),
        }
    )

    global_score = float(sum(float(c["score"]) for c in checks) / len(checks))

    if pta_ready and transfer_pass and not gw_redesign_required and not structural_redesign:
        readiness = "GLOBAL_EMPIRICAL_READY"
        hard_gate = "PASS"
        recommendation = "START_JOINT_PTA_GW_CONFIRMATORY_CAMPAIGN"
    elif pta_ready and transfer_partial and (gw_redesign_required or structural_redesign):
        readiness = "GLOBAL_PARTIAL_PTA_ONLY"
        hard_gate = "PARTIAL"
        recommendation = "RUN_PTA_TRACK_CONTINUATION_AND_PARALLEL_GW_REDESIGN_PROGRAM"
    else:
        readiness = "GLOBAL_NOT_READY"
        hard_gate = "FAIL"
        recommendation = "HOLD_GLOBAL_CLAIMS_AND_FOCUS_ON_METHOD_REPAIR"

    next_steps = [
        "QW-1831: event-windowed GW coherent feature extractor (chirp-conditioned).",
        "QW-1832: GW shared-vs-control objective with explicit near-target prevalence optimization.",
        "QW-1833: GW multi-detector consistency gate on H1-L1, H1-V1, L1-V1 with quantile scoring.",
    ]

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "global_score": global_score,
        "hard_gate": hard_gate,
        "readiness": readiness,
        "recommendation": recommendation,
        "next_steps": next_steps,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1830: GLOBAL EMPIRICAL STATUS GATE",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Global score: {global_score:.3f}",
        f"- Hard gate: **{hard_gate}**",
        f"- Readiness: **{readiness}**",
        f"- Recommendation: **{recommendation}**",
        "",
        "## Checks",
    ]
    for c in checks:
        lines.append(f"- {c['domain']}: {c['status']} | score={c['score']:.3f} | note={c['note']}")

    lines.extend(["", "## Next Steps"])
    for s in next_steps:
        lines.append(f"- {s}")

    lines.extend(["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1830] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1830] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
