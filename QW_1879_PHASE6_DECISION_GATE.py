#!/usr/bin/env python3
"""
QW-1879: Phase-6 decision gate.

Final decision after phase-5 and reformulation branch.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1879_phase6_decision_gate.json"
OUT_MD = ROOT / "RAPORT_QW1879_PHASE6_DECISION_GATE.md"


def read_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def clip01(x: float) -> float:
    return max(0.0, min(1.0, float(x)))


def main() -> None:
    d1876 = read_json("report_qw1876_kernel_closure_phase5_gate.json")
    d1878 = read_json("report_qw1878_reformulation_mechanism_comparison.json")

    s5 = clip01(float(d1876.get("global_score", 0.0)))
    s_cmp = clip01(float(d1878.get("comparison_score", 0.0)))

    decision_score = 0.60 * s5 + 0.40 * s_cmp

    checks = {
        "phase5_closure": {
            "score": s5,
            "pass": s5 >= 0.62,
            "note": d1876.get("readiness"),
        },
        "reformulation_superiority": {
            "score": s_cmp,
            "pass": s_cmp >= 0.65,
            "note": d1878.get("verdict"),
        },
    }

    if all(v["pass"] for v in checks.values()) and decision_score >= 0.68:
        readiness = "INTERNAL_CLOSURE_PATH_AVAILABLE"
    elif decision_score >= 0.50:
        readiness = "INTERNAL_PATH_PARTIAL_REQUIRES_NEW_HYPOTHESIS"
    else:
        readiness = "INTERNAL_CLOSURE_NOT_AVAILABLE_UNDER_CURRENT_MICROMODEL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "decision_score": decision_score,
        "readiness": readiness,
        "verdict": "PHASE6_DECISION_GATE_COMPLETE",
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1879: PHASE-6 DECISION GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        f"- Decision score: {decision_score:.3f}",
        f"- Readiness: **{readiness}**",
        "",
        "## Checks",
    ]

    for k, v in checks.items():
        lines.append(f"- {k}: {'PASS' if v['pass'] else 'FAIL'} | score={v['score']:.3f} | note={v['note']}")

    lines += ["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1879] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1879] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
