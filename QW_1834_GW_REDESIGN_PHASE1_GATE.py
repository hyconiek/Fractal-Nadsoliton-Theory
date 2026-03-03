#!/usr/bin/env python3
"""
QW-1834: GW redesign phase-1 progress gate.

Aggregates legacy failures and new redesign outcomes:
- QW-1827: GW redesign gate (legacy diagnostics),
- QW-1829: near-target feasibility,
- QW-1832: shared-control window objective,
- QW-1833: multi-detector consistency gate.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1834_gw_redesign_phase1_gate.json"
OUT_MD = ROOT / "RAPORT_QW1834_GW_REDESIGN_PHASE1_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    d1827 = load("report_qw1827_gw_redesign_gate.json")
    d1829 = load("report_qw1829_gw_near_target_feasibility.json")
    d1832 = load("report_qw1832_gw_shared_control_window_objective.json")
    d1833 = load("report_qw1833_gw_multi_detector_consistency_gate.json")

    checks: List[Dict[str, object]] = []

    # Legacy baseline severity (the lower 1827 score, the lower this contribution)
    s_legacy = float(d1827.get("global_score", 0.0))
    checks.append(
        {
            "domain": "Legacy GW baseline status (1827)",
            "score": s_legacy,
            "status": "PASS" if s_legacy >= 0.60 else "FAIL",
            "note": d1827.get("readiness", ""),
        }
    )

    # Structural feasibility remains a hard blocker if verdict is structural redesign
    hard_struct = d1829.get("verdict") == "GW_NEAR_TARGET_REQUIRES_STRUCTURAL_REDESIGN"
    s_feas = 0.0 if hard_struct else 0.8
    checks.append(
        {
            "domain": "Near-target structural feasibility (1829)",
            "score": s_feas,
            "status": "FAIL" if hard_struct else "PASS",
            "note": d1829.get("verdict", ""),
        }
    )

    p1832 = d1832.get("pass_flags", {})
    s_obj = (
        0.40 * float(bool(p1832.get("auc_support")))
        + 0.30 * float(bool(p1832.get("balanced_accuracy_support")))
        + 0.30 * float(bool(p1832.get("prevalence_advantage_support")))
    )
    checks.append(
        {
            "domain": "Shared-control objective support (1832)",
            "score": s_obj,
            "status": "PASS" if s_obj >= 0.90 else "FAIL",
            "note": d1832.get("verdict", ""),
        }
    )

    p1833 = d1833.get("pass_flags", {})
    s_cons = (
        0.25 * float(bool(p1833.get("shared_vs_ctrl_separation")))
        + 0.25 * float(bool(p1833.get("shared_prevalence_advantage")))
        + 0.25 * float(bool(p1833.get("auc_support")))
        + 0.25 * float(bool(p1833.get("control_pair_consistency")))
    )
    checks.append(
        {
            "domain": "Multi-detector consistency (1833)",
            "score": s_cons,
            "status": "PASS" if s_cons >= 0.75 else "FAIL",
            "note": d1833.get("verdict", ""),
        }
    )

    global_score = float(sum(float(c["score"]) for c in checks) / len(checks))

    # Progress logic: require both new modules strong and no hard structural blocker for readiness.
    new_modules_strong = (s_obj >= 0.90) and (s_cons >= 0.75)

    if new_modules_strong and not hard_struct:
        readiness = "GW_REDESIGN_PHASE1_READY"
        hard_gate = "PASS"
        recommendation = "MOVE_TO_GW_CONFIRMATORY_PHASE"
    elif new_modules_strong and hard_struct:
        readiness = "GW_REDESIGN_PHASE1_PARTIAL_PROGRESS"
        hard_gate = "PARTIAL"
        recommendation = "KEEP_NEW_OBJECTIVE_AND_TARGET_STRUCTURAL_REPARAM_NEXT"
    else:
        readiness = "GW_REDESIGN_PHASE1_NOT_READY"
        hard_gate = "FAIL"
        recommendation = "CONTINUE_GW_METHOD_REBUILD"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "global_score": global_score,
        "hard_gate": hard_gate,
        "readiness": readiness,
        "recommendation": recommendation,
        "new_module_scores": {
            "objective_1832": s_obj,
            "consistency_1833": s_cons,
        },
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1834: GW REDESIGN PHASE-1 GATE",
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

    lines.extend(
        [
            "",
            "## New Module Scores",
            f"- objective_1832: {s_obj:.3f}",
            f"- consistency_1833: {s_cons:.3f}",
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1834] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1834] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
