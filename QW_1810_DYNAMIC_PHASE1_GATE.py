#!/usr/bin/env python3
"""
QW-1810: Dynamic phase-1 gate.

Aggregates first dynamic attempts (QW-1808, QW-1809) after transition gate (QW-1807).
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1810_dynamic_phase1_gate.json"
OUT_MD = ROOT / "RAPORT_QW1810_DYNAMIC_PHASE1_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    d1807 = load("report_qw1807_next_phase_gate.json")
    d1808 = load("report_qw1808_dynamic_drift_regime_model.json")
    d1809 = load("report_qw1809_dynamic_feature_scan.json")

    checks: List[Dict[str, object]] = []

    # Phase transition precondition
    pass07 = d1807.get("readiness") == "TRANSITION_TO_DYNAMIC_PHASE_RECOMMENDED"
    checks.append(
        {
            "domain": "Transition gate precondition (1807)",
            "score": 1.0 if pass07 else 0.0,
            "status": "PASS" if pass07 else "FAIL",
            "note": d1807.get("recommendation", ""),
        }
    )

    p08 = d1808.get("pass_flags", {})
    s08 = (
        (0.30 if p08.get("full_gain") else 0.0)
        + (0.30 if p08.get("replication_gain") else 0.0)
        + (0.20 if p08.get("dispersion_control") else 0.0)
        + (0.20 if p08.get("dynamic_residual_reduction") else 0.0)
    )
    checks.append(
        {
            "domain": "Single dynamic drift model (1808)",
            "score": s08,
            "status": "PASS" if s08 >= 0.75 else "FAIL",
            "note": d1808.get("verdict", ""),
        }
    )

    p09 = d1809.get("pass_flags", {})
    s09 = (
        (0.55 if p09.get("gain") else 0.0)
        + (0.45 if p09.get("dispersion_control") else 0.0)
    )
    checks.append(
        {
            "domain": "Dynamic feature scan (1809)",
            "score": s09,
            "status": "PASS" if s09 >= 0.75 else "FAIL",
            "note": d1809.get("verdict", ""),
        }
    )

    global_score = float(sum(float(c["score"]) for c in checks) / len(checks))
    hard_gate = all(c["status"] == "PASS" for c in checks)

    if hard_gate:
        readiness = "DYNAMIC_PHASE1_READY"
        recommendation = "CONTINUE_SIMPLE_DYNAMIC_MODELS"
    elif global_score >= 0.45:
        readiness = "DYNAMIC_PHASE1_PARTIAL"
        recommendation = "ESCALATE_TO_DYNAMIC_PHASE2_LATENT_STATE"
    else:
        readiness = "DYNAMIC_PHASE1_NOT_READY"
        recommendation = "REDESIGN_DYNAMIC_INPUT_REPRESENTATION"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "global_score": global_score,
        "hard_gate": "PASS" if hard_gate else "FAIL",
        "readiness": readiness,
        "recommendation": recommendation,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1810: DYNAMIC PHASE-1 GATE",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Global score: {global_score:.3f}",
        f"- Hard gate: **{output['hard_gate']}**",
        f"- Readiness: **{readiness}**",
        f"- Recommendation: **{recommendation}**",
        "",
        "## Checks",
    ]
    for c in checks:
        lines.append(f"- {c['domain']}: {c['status']} | score={c['score']:.3f} | note={c['note']}")
    lines.extend(["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1810] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1810] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
