#!/usr/bin/env python3
"""
QW-1800: Hierarchical branch gate.

Aggregates QW-1796..1799 to decide whether hierarchical multimode branch
is operationally ready or should be parked.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1800_hierarchical_branch_gate.json"
OUT_MD = ROOT / "RAPORT_QW1800_HIERARCHICAL_BRANCH_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    d1796 = load("report_qw1796_physical_multimode_extension.json")
    d1797 = load("report_qw1797_hierarchical_multimode_shrinkage.json")
    d1798 = load("report_qw1798_hierarchical_shrinkage_calibration.json")
    d1799 = load("report_qw1799_hierarchical_transfer_test.json")

    checks: List[Dict[str, object]] = []

    # 1796: branch seed quality.
    p96 = d1796.get("pass_flags", {})
    s96 = (
        (0.30 if p96.get("full_gain_over_m2") else 0.0)
        + (0.30 if p96.get("replication_gain") else 0.0)
        + (0.20 if p96.get("dispersion_control") else 0.0)
        + (0.20 if p96.get("residual_flattening") else 0.0)
    )
    checks.append(
        {
            "domain": "Seed multimode extension (1796)",
            "score": s96,
            "status": "PASS" if s96 >= 0.75 else "FAIL",
            "note": d1796.get("verdict", ""),
        }
    )

    # 1797: hierarchical gain.
    p97 = d1797.get("pass_flags", {})
    s97 = (
        (0.30 if p97.get("full_gain_over_m2") else 0.0)
        + (0.30 if p97.get("replication_gain") else 0.0)
        + (0.20 if p97.get("dispersion_control") else 0.0)
        + (0.20 if p97.get("residual_improvement") else 0.0)
    )
    checks.append(
        {
            "domain": "Hierarchical shrinkage gain (1797)",
            "score": s97,
            "status": "PASS" if s97 >= 0.75 else "FAIL",
            "note": d1797.get("verdict", ""),
        }
    )

    # 1798: calibration quality.
    best98 = d1798.get("best", {})
    p98 = bool(best98.get("pass_basic", False))
    s98 = 1.0 if p98 else 0.55 if d1798.get("recommendation_strength") == "CONDITIONAL" else 0.0
    checks.append(
        {
            "domain": "Shrinkage calibration (1798)",
            "score": s98,
            "status": "PASS" if s98 >= 0.75 else "FAIL",
            "note": d1798.get("verdict", ""),
        }
    )

    # 1799: transfer test is decisive.
    p99 = d1799.get("pass_flags", {})
    s99 = (
        (0.35 if p99.get("transfer_gain") else 0.0)
        + (0.25 if p99.get("transfer_stability") else 0.0)
        + (0.20 if p99.get("generalization_gap_control") else 0.0)
        + (0.20 if p99.get("test_positive_vs_flat") else 0.0)
    )
    checks.append(
        {
            "domain": "Transfer robustness (1799)",
            "score": s99,
            "status": "PASS" if s99 >= 0.80 else "FAIL",
            "note": d1799.get("verdict", ""),
        }
    )

    global_score = float(sum(float(c["score"]) for c in checks) / len(checks))
    hard_gate = all(c["status"] == "PASS" for c in checks)

    if hard_gate and global_score >= 0.85:
        readiness = "HIERARCHICAL_BRANCH_READY"
    elif global_score >= 0.60:
        readiness = "HIERARCHICAL_BRANCH_PARTIAL"
    else:
        readiness = "HIERARCHICAL_BRANCH_NOT_READY"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "global_score": global_score,
        "hard_gate": "PASS" if hard_gate else "FAIL",
        "readiness": readiness,
        "recommendation": (
            "CONTINUE_WITH_HIERARCHICAL_BRANCH" if hard_gate else "PARK_BRANCH_AND_REQUIRE_NEW_PHYSICAL_MECHANISM"
        ),
    }

    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1800: HIERARCHICAL BRANCH GATE",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Global score: {global_score:.3f}",
        f"- Hard gate: **{output['hard_gate']}**",
        f"- Readiness: **{readiness}**",
        f"- Recommendation: **{output['recommendation']}**",
        "",
        "## Checks",
    ]
    for c in checks:
        lines.append(f"- {c['domain']}: {c['status']} | score={c['score']:.3f} | note={c['note']}")
    lines.extend(["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1800] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1800] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
