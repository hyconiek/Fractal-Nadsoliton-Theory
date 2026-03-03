#!/usr/bin/env python3
"""
QW-1814: Dynamic phase-2 gate.

Aggregates dynamic gates/results:
- QW-1810 (dynamic phase-1 gate),
- QW-1812 (dynamic proxy branch gate),
- QW-1813 (sequence-window model).
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1814_dynamic_phase2_gate.json"
OUT_MD = ROOT / "RAPORT_QW1814_DYNAMIC_PHASE2_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    d1810 = load("report_qw1810_dynamic_phase1_gate.json")
    d1812 = load("report_qw1812_dynamic_proxy_branch_gate.json")
    d1813 = load("report_qw1813_sequence_window_model.json")

    checks: List[Dict[str, object]] = []

    s10 = float(d1810.get("global_score", 0.0))
    pass10 = d1810.get("hard_gate") == "PASS"
    checks.append(
        {
            "domain": "Dynamic phase-1 gate (1810)",
            "score": s10,
            "status": "PASS" if pass10 else "FAIL",
            "note": d1810.get("readiness", ""),
        }
    )

    s12 = float(d1812.get("global_score", 0.0))
    pass12 = d1812.get("hard_gate") == "PASS"
    checks.append(
        {
            "domain": "Dynamic proxy branch gate (1812)",
            "score": s12,
            "status": "PASS" if pass12 else "FAIL",
            "note": d1812.get("readiness", ""),
        }
    )

    p13 = d1813.get("pass_flags", {})
    s13 = (
        (0.40 if p13.get("full_gain") else 0.0)
        + (0.35 if p13.get("replication_gain") else 0.0)
        + (0.25 if p13.get("dispersion_control") else 0.0)
    )
    checks.append(
        {
            "domain": "Sequence-window model (1813)",
            "score": s13,
            "status": "PASS" if s13 >= 0.75 else "FAIL",
            "note": d1813.get("verdict", ""),
        }
    )

    global_score = float(sum(float(c["score"]) for c in checks) / len(checks))
    hard_gate = all(c["status"] == "PASS" for c in checks)

    if hard_gate:
        readiness = "DYNAMIC_PHASE2_READY"
        recommendation = "CONTINUE_DYNAMIC_MODEL_DEVELOPMENT"
    elif global_score >= 0.45:
        readiness = "DYNAMIC_PHASE2_PARTIAL"
        recommendation = "COLLECT_RICHER_SEQUENCE_DATA_BEFORE_NEXT_MODEL_ITERATION"
    else:
        readiness = "DYNAMIC_PHASE2_NOT_READY"
        recommendation = "PARK_DYNAMIC_BRANCH_AND_REDESIGN_DATA_REPRESENTATION_PIPELINE"

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
        "# RAPORT QW-1814: DYNAMIC PHASE-2 GATE",
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

    print(f"[QW-1814] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1814] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
