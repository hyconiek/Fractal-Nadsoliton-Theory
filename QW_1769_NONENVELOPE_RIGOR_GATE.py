#!/usr/bin/env python3
"""
QW-1769: Non-envelope rigor gate (1764, 1766, 1768).
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1769_nonenvelope_rigor_gate.json"
OUT_MD = ROOT / "RAPORT_QW1769_NONENVELOPE_RIGOR_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    d1764 = load("report_qw1764_memory_field_gate.json")
    d1766 = load("report_qw1766_nonenvelope_paired_intervention.json")
    d1768 = load("report_qw1768_leakage_controlled_paired_intervention.json")

    checks: List[Dict[str, object]] = []

    s64 = float(d1764.get("global_score", 0.0))
    checks.append(
        {
            "domain": "Legacy envelope branch baseline (1764)",
            "score": s64,
            "status": "PASS" if s64 >= 0.60 else "FAIL",
            "note": d1764.get("readiness", ""),
        }
    )

    p66 = d1766.get("pass_flags", {})
    s66 = (
        (0.30 if p66.get("discrimination") else 0.0)
        + (0.30 if p66.get("continuous_regression") else 0.0)
        + (0.20 if p66.get("calibration_nonboundary") else 0.0)
        + (0.20 if p66.get("cv_stability") else 0.0)
    )
    checks.append(
        {
            "domain": "Paired intervention (non-grouped CV) 1766",
            "score": s66,
            "status": "PASS" if s66 >= 0.70 else "FAIL",
            "note": d1766.get("verdict", ""),
        }
    )

    p68 = d1768.get("pass_flags", {})
    s68 = (
        (0.30 if p68.get("discrimination") else 0.0)
        + (0.30 if p68.get("continuous_regression") else 0.0)
        + (0.20 if p68.get("calibration_nonboundary") else 0.0)
        + (0.20 if p68.get("cv_stability") else 0.0)
    )
    checks.append(
        {
            "domain": "Leakage-controlled paired intervention 1768",
            "score": s68,
            "status": "PASS" if s68 >= 0.70 else "FAIL",
            "note": d1768.get("verdict", ""),
        }
    )

    global_score = float(sum(float(c["score"]) for c in checks) / len(checks))
    hard_gate = all(c["status"] == "PASS" for c in checks)

    if hard_gate and global_score >= 0.80:
        readiness = "NONENVELOPE_RIGOR_CLOSED"
    elif global_score >= 0.60:
        readiness = "NONENVELOPE_RIGOR_PARTIAL"
    else:
        readiness = "NONENVELOPE_RIGOR_OPEN"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "global_score": global_score,
        "hard_gate": "PASS" if hard_gate else "FAIL",
        "readiness": readiness,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1769: NONENVELOPE RIGOR GATE",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Global score: {global_score:.3f}",
        f"- Hard gate: **{output['hard_gate']}**",
        f"- Readiness: **{readiness}**",
        "",
        "## Checks",
    ]
    for c in checks:
        lines.append(f"- {c['domain']}: {c['status']} | score={c['score']:.3f} | note={c['note']}")
    lines.extend(["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1769] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1769] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
