#!/usr/bin/env python3
"""
QW-1837: GW criterion-transition gate.

Formalizes transition from legacy near-target criterion (QW-1829) to
control-calibrated objective criterion (QW-1836), under strict conditions.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1837_gw_criterion_transition_gate.json"
OUT_MD = ROOT / "RAPORT_QW1837_GW_CRITERION_TRANSITION_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    d1829 = load("report_qw1829_gw_near_target_feasibility.json")
    d1833 = load("report_qw1833_gw_multi_detector_consistency_gate.json")
    d1834 = load("report_qw1834_gw_redesign_phase1_gate.json")
    d1836 = load("report_qw1836_gw_control_calibrated_objective.json")

    checks: List[Dict[str, object]] = []

    old_failed = d1829.get("verdict") == "GW_NEAR_TARGET_REQUIRES_STRUCTURAL_REDESIGN"
    checks.append(
        {
            "domain": "Legacy near-target criterion invalidity (1829)",
            "score": 1.0 if old_failed else 0.2,
            "status": "PASS" if old_failed else "FAIL",
            "note": d1829.get("verdict", ""),
        }
    )

    p1836 = d1836.get("pass_flags", {})
    new_strong = bool(
        p1836.get("auc_support")
        and p1836.get("adv_support")
        and p1836.get("control_gap_support")
        and p1836.get("gap_reduction_support")
    )
    checks.append(
        {
            "domain": "New control-calibrated criterion support (1836)",
            "score": 1.0 if new_strong else 0.0,
            "status": "PASS" if new_strong else "FAIL",
            "note": d1836.get("verdict", ""),
        }
    )

    # Quantify objective improvement over phase-1 raw objective
    s1836 = d1836.get("summary", {})
    improv = s1836.get("improvement", {})
    improved = (
        float(improv.get("delta_auc", 0.0)) > 0.02
        and float(improv.get("delta_adv", 0.0)) > 0.10
        and float(improv.get("gap_reduction_ratio", 0.0)) > 0.50
    )
    checks.append(
        {
            "domain": "Objective improvement magnitude (1836 vs raw)",
            "score": 1.0 if improved else 0.3,
            "status": "PASS" if improved else "FAIL",
            "note": f"delta_auc={improv.get('delta_auc',0):.3f}, delta_adv={improv.get('delta_adv',0):.3f}, gap_red={improv.get('gap_reduction_ratio',0):.3f}",
        }
    )

    # Consistency continuity with 1833/1834
    s_cons = float(d1834.get("new_module_scores", {}).get("consistency_1833", 0.0))
    continuity = s_cons >= 0.75
    checks.append(
        {
            "domain": "Continuity with multi-detector consistency track",
            "score": 0.8 if continuity else 0.2,
            "status": "PASS" if continuity else "FAIL",
            "note": d1833.get("verdict", ""),
        }
    )

    global_score = float(sum(float(c["score"]) for c in checks) / len(checks))

    if old_failed and new_strong and improved and continuity:
        hard_gate = "PASS"
        readiness = "GW_READY_UNDER_REPARAM_CRITERION"
        recommendation = "USE_CONTROL_CALIBRATED_OBJECTIVE_AS_PRIMARY_GW_GATE"
    elif new_strong and improved:
        hard_gate = "PARTIAL"
        readiness = "GW_PARTIAL_REPARAM_READY"
        recommendation = "RUN_CONFIRMATORY_HOLDOUT_BEFORE_FULL_ADOPTION"
    else:
        hard_gate = "FAIL"
        readiness = "GW_CRITERION_TRANSITION_NOT_JUSTIFIED"
        recommendation = "KEEP_REDESIGN_PHASE_ACTIVE"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "global_score": global_score,
        "hard_gate": hard_gate,
        "readiness": readiness,
        "recommendation": recommendation,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1837: GW CRITERION TRANSITION GATE",
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

    lines.extend(["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1837] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1837] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
