#!/usr/bin/env python3
"""
QW-1745: Aggregate rigor gate for micromodel iteration (1739-1744).
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1745_micromodel_rigor_gate.json"
OUT_MD = ROOT / "RAPORT_QW1745_MICROMODEL_RIGOR_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    d1739 = load("report_qw1739_signed_dynamic_micromodel_derivation.json")
    d1740 = load("report_qw1740_signed_micromodel_identifiability_audit.json")
    d1741 = load("report_qw1741_constrained_global_derivation.json")
    d1742 = load("report_qw1742_profile_likelihood_identifiability.json")
    d1743 = load("report_qw1743_oscillatory_cohort_derivation.json")
    d1744 = load("report_qw1744_oscillatory_cohort_identifiability.json")

    checks: List[Dict[str, object]] = []

    # 1739
    p1739 = d1739["pass_flags"]
    s1739 = (0.5 if p1739.get("fit_quality") else 0.0) + (0.3 if p1739.get("stability") else 0.0) + (0.2 if p1739.get("near_reference") else 0.0)
    checks.append(
        {
            "domain": "Signed dynamic derivation (1739)",
            "score": s1739,
            "status": "PASS" if s1739 >= 0.7 else "FAIL",
            "note": d1739["verdict"],
        }
    )

    # 1740
    p1740 = d1740["pass_flags"]
    s1740 = (0.35 if p1740.get("conditioning") else 0.0) + (0.30 if p1740.get("ci_width") else 0.0) + (0.20 if p1740.get("mode_count") else 0.0) + (0.15 if p1740.get("correlation") else 0.0)
    checks.append(
        {
            "domain": "Identifiability audit (1740)",
            "score": s1740,
            "status": "PASS" if s1740 >= 0.7 else "FAIL",
            "note": d1740["verdict"],
        }
    )

    # 1741
    p1741 = d1741["pass_flags"]
    s1741 = (0.45 if p1741.get("fit_quality") else 0.0) + (0.30 if p1741.get("identifiability") else 0.0) + (0.25 if p1741.get("nonboundary_solution") else 0.0)
    checks.append(
        {
            "domain": "Constrained global fit (1741)",
            "score": s1741,
            "status": "PASS" if s1741 >= 0.7 else "FAIL",
            "note": d1741["verdict"],
        }
    )

    # 1742
    p1742 = d1742["pass_flags"]
    s1742 = (0.45 if p1742.get("profile_width") else 0.0) + (0.30 if p1742.get("conditioning") else 0.0) + (0.25 if p1742.get("correlation") else 0.0)
    checks.append(
        {
            "domain": "Profile-likelihood identifiability (1742)",
            "score": s1742,
            "status": "PASS" if s1742 >= 0.7 else "FAIL",
            "note": d1742["verdict"],
        }
    )

    # 1743
    p1743 = d1743["pass_flags"]
    s1743 = (0.35 if p1743.get("fit_quality") else 0.0) + (0.25 if p1743.get("informative_cohort") else 0.0) + (0.20 if p1743.get("nonboundary_solution") else 0.0) + (0.20 if p1743.get("ci_width") else 0.0)
    checks.append(
        {
            "domain": "Oscillatory cohort derivation (1743)",
            "score": s1743,
            "status": "PASS" if s1743 >= 0.7 else "FAIL",
            "note": d1743["verdict"],
        }
    )

    # 1744
    p1744 = d1744["pass_flags"]
    s1744 = (0.45 if p1744.get("profile_width") else 0.0) + (0.30 if p1744.get("conditioning") else 0.0) + (0.25 if p1744.get("correlation") else 0.0)
    checks.append(
        {
            "domain": "Oscillatory cohort identifiability (1744)",
            "score": s1744,
            "status": "PASS" if s1744 >= 0.7 else "FAIL",
            "note": d1744["verdict"],
        }
    )

    global_score = float(sum(float(c["score"]) for c in checks) / len(checks))
    hard_gate = all(c["status"] == "PASS" for c in checks)

    if hard_gate and global_score >= 0.8:
        readiness = "MICROMODEL_ITERATION_CLOSED"
    elif global_score >= 0.6:
        readiness = "MICROMODEL_ITERATION_PARTIAL"
    else:
        readiness = "MICROMODEL_ITERATION_OPEN"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "global_score": global_score,
        "hard_gate": "PASS" if hard_gate else "FAIL",
        "readiness": readiness,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1745: MICROMODEL RIGOR GATE",
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

    print(f"[QW-1745] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1745] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
