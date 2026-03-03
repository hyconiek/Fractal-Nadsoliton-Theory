#!/usr/bin/env python3
"""
QW-1838: Global readiness gate under reparameterized GW criterion.

Updates global status by replacing legacy GW near-target criterion with
QW-1837-approved control-calibrated GW criterion.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1838_global_reparam_readiness_gate.json"
OUT_MD = ROOT / "RAPORT_QW1838_GLOBAL_REPARAM_READINESS_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    d1824 = load("report_qw1824_quantile_gated_readiness.json")
    d1825 = load("report_qw1825_quantile_protocol_transfer_gate.json")
    d1836 = load("report_qw1836_gw_control_calibrated_objective.json")
    d1837 = load("report_qw1837_gw_criterion_transition_gate.json")
    d1835 = load("report_qw1835_global_status_after_gw_phase1.json")

    checks: List[Dict[str, object]] = []

    pta_ready = d1824.get("hard_gate") == "PASS"
    checks.append(
        {
            "domain": "PTA branch readiness (1824)",
            "score": 1.0 if pta_ready else 0.0,
            "status": "PASS" if pta_ready else "FAIL",
            "note": d1824.get("readiness", ""),
        }
    )

    gw_objective_ready = d1836.get("verdict") == "GW_CONTROL_CALIBRATED_OBJECTIVE_SUPPORTED"
    checks.append(
        {
            "domain": "GW control-calibrated objective (1836)",
            "score": 1.0 if gw_objective_ready else 0.0,
            "status": "PASS" if gw_objective_ready else "FAIL",
            "note": d1836.get("verdict", ""),
        }
    )

    gw_transition_ready = d1837.get("hard_gate") == "PASS"
    checks.append(
        {
            "domain": "GW criterion-transition justification (1837)",
            "score": 1.0 if gw_transition_ready else 0.0,
            "status": "PASS" if gw_transition_ready else "FAIL",
            "note": d1837.get("readiness", ""),
        }
    )

    legacy_transfer = d1825.get("hard_gate", "FAIL")
    legacy_score = 1.0 if legacy_transfer == "PASS" else (0.6 if legacy_transfer == "PARTIAL" else 0.0)
    checks.append(
        {
            "domain": "Legacy transfer status (1825) [for traceability]",
            "score": legacy_score,
            "status": legacy_transfer,
            "note": d1825.get("readiness", ""),
        }
    )

    old_global = float(d1835.get("global_score_updated", 0.0))
    new_score = float(sum(float(c["score"]) for c in checks) / len(checks))
    delta = new_score - old_global

    if pta_ready and gw_objective_ready and gw_transition_ready:
        hard_gate = "PASS"
        readiness = "GLOBAL_CONDITIONAL_READY_UNDER_REPARAM_CRITERIA"
        recommendation = "START_PRE_REGISTERED_JOINT_CONFIRMATORY_CAMPAIGN_PTA_GW"
    elif pta_ready and (gw_objective_ready or gw_transition_ready):
        hard_gate = "PARTIAL"
        readiness = "GLOBAL_PARTIAL_REPARAM"
        recommendation = "COMPLETE_GW_REPARAM_VALIDATION_AND_RETEST"
    else:
        hard_gate = "FAIL"
        readiness = "GLOBAL_NOT_READY"
        recommendation = "CONTINUE_METHOD_REPAIR"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "global_score_previous_qw1835": old_global,
        "global_score_updated": new_score,
        "delta_global_score": delta,
        "hard_gate": hard_gate,
        "readiness": readiness,
        "recommendation": recommendation,
        "note": "Readiness is conditional on reparameterized GW criterion replacing legacy near-target criterion.",
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1838: GLOBAL REPARAM READINESS GATE",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Global score (QW-1835 -> now): {old_global:.3f} -> {new_score:.3f}",
        f"- Delta score: {delta:.3f}",
        f"- Hard gate: **{hard_gate}**",
        f"- Readiness: **{readiness}**",
        f"- Recommendation: **{recommendation}**",
        f"- Note: {output['note']}",
        "",
        "## Checks",
    ]
    for c in checks:
        lines.append(f"- {c['domain']}: {c['status']} | score={c['score']:.3f} | note={c['note']}")

    lines.extend(["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1838] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1838] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
