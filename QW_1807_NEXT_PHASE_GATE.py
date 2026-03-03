#!/usr/bin/env python3
"""
QW-1807: Next-phase gate after static-branch exhaustion.

Purpose:
- consolidate branch-level gates and diagnostics,
- decide whether to continue static angular/model-class tuning
  or transition to dynamic latent-regime modeling.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1807_next_phase_gate.json"
OUT_MD = ROOT / "RAPORT_QW1807_NEXT_PHASE_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    d1793 = load("report_qw1793_model_lockin_gate.json")
    d1800 = load("report_qw1800_hierarchical_branch_gate.json")
    d1805 = load("report_qw1805_decoherence_branch_gate.json")
    d1806 = load("report_qw1806_quality_regime_diagnostic.json")

    checks: List[Dict[str, object]] = []

    # Baseline lock-in must stay valid.
    s93 = float(d1793.get("global_score", 0.0))
    pass93 = d1793.get("hard_gate") == "PASS" and s93 >= 0.90
    checks.append(
        {
            "domain": "Locked baseline validity (1793)",
            "score": s93,
            "status": "PASS" if pass93 else "FAIL",
            "note": d1793.get("readiness", ""),
        }
    )

    # Hierarchical branch readiness.
    s00 = float(d1800.get("global_score", 0.0))
    pass00 = d1800.get("hard_gate") == "PASS"
    checks.append(
        {
            "domain": "Hierarchical static branch (1800)",
            "score": s00,
            "status": "PASS" if pass00 else "FAIL",
            "note": d1800.get("readiness", ""),
        }
    )

    # Decoherence branch readiness.
    s05 = float(d1805.get("global_score", 0.0))
    pass05 = d1805.get("hard_gate") == "PASS"
    checks.append(
        {
            "domain": "Decoherence static branch (1805)",
            "score": s05,
            "status": "PASS" if pass05 else "FAIL",
            "note": d1805.get("readiness", ""),
        }
    )

    # Quality-regime diagnostic should not support static quality correction as main mechanism.
    s06 = d1806.get("summary", {})
    delta06 = float(s06.get("full_delta_m2q_vs_m2", 0.0))
    pass06 = delta06 <= 0.0 and d1806.get("verdict") == "QUALITY_REGIME_SIGNAL_WEAK"
    checks.append(
        {
            "domain": "Quality-regime correction check (1806)",
            "score": 1.0 if pass06 else 0.0,
            "status": "PASS" if pass06 else "FAIL",
            "note": d1806.get("verdict", ""),
        }
    )

    global_score = float(sum(float(c["score"]) for c in checks) / len(checks))
    static_ready = pass00 or pass05

    if pass93 and (not static_ready) and pass06:
        readiness = "TRANSITION_TO_DYNAMIC_PHASE_RECOMMENDED"
        recommendation = "LAUNCH_DYNAMIC_LATENT_REGIME_PROGRAM"
    elif static_ready:
        readiness = "STATIC_PHASE_CONTINUE_OPTIONAL"
        recommendation = "CONTINUE_BEST_STATIC_BRANCH"
    else:
        readiness = "INSUFFICIENT_SIGNAL_FOR_TRANSITION"
        recommendation = "REVIEW_ASSUMPTIONS"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "global_score": global_score,
        "static_branch_ready": bool(static_ready),
        "readiness": readiness,
        "recommendation": recommendation,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1807: NEXT PHASE GATE",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Global score: {global_score:.3f}",
        f"- Static branch ready: {output['static_branch_ready']}",
        f"- Readiness: **{readiness}**",
        f"- Recommendation: **{recommendation}**",
        "",
        "## Checks",
    ]
    for c in checks:
        lines.append(f"- {c['domain']}: {c['status']} | score={c['score']:.3f} | note={c['note']}")
    lines.extend(["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1807] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1807] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
