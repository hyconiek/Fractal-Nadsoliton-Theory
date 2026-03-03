#!/usr/bin/env python3
"""
QW-1805: Decoherence branch gate.

Aggregates QW-1801..1804 to decide operational status of the
physical decoherence extension branch.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1805_decoherence_branch_gate.json"
OUT_MD = ROOT / "RAPORT_QW1805_DECOHERENCE_BRANCH_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    d1801 = load("report_qw1801_physical_decoherence_extension.json")
    d1802 = load("report_qw1802_decoherence_form_comparison.json")
    d1803 = load("report_qw1803_additive_hierarchical_stabilization.json")
    d1804 = load("report_qw1804_variance_decomposition_mc.json")

    checks: List[Dict[str, object]] = []

    p01 = d1801.get("pass_flags", {})
    s01 = (
        (0.25 if p01.get("full_gain_over_m2") else 0.0)
        + (0.25 if p01.get("replication_gain") else 0.0)
        + (0.25 if p01.get("dispersion_control") else 0.0)
        + (0.25 if p01.get("residual_flattening") else 0.0)
    )
    checks.append(
        {
            "domain": "Seed decoherence extension (1801)",
            "score": s01,
            "status": "PASS" if s01 >= 0.75 else "FAIL",
            "note": d1801.get("verdict", ""),
        }
    )

    p02 = d1802.get("pass_flags", {})
    s02 = (
        (0.30 if p02.get("additive_gain_over_m2") else 0.0)
        + (0.20 if p02.get("additive_better_than_multiplicative") else 0.0)
        + (0.25 if p02.get("additive_dispersion_control") else 0.0)
        + (0.25 if p02.get("additive_residual_flattening") else 0.0)
    )
    checks.append(
        {
            "domain": "Form comparison (1802)",
            "score": s02,
            "status": "PASS" if s02 >= 0.75 else "FAIL",
            "note": d1802.get("verdict", ""),
        }
    )

    best03 = d1803.get("best", {})
    pass03 = bool(best03.get("pass_basic", False))
    s03 = 1.0 if pass03 else 0.60 if d1803.get("recommendation_strength") == "CONDITIONAL" else 0.0
    checks.append(
        {
            "domain": "Additive hierarchical stabilization (1803)",
            "score": s03,
            "status": "PASS" if s03 >= 0.75 else "FAIL",
            "note": d1803.get("verdict", ""),
        }
    )

    p04 = d1804.get("pass_flags", {})
    # Here we *require* manageable split stability; MC-dominance is optional.
    s04 = (
        (0.30 if p04.get("positive_mean_delta") else 0.0)
        + (0.20 if p04.get("mc_dominant_component") else 0.0)
        + (0.50 if p04.get("between_split_std_small") else 0.0)
    )
    checks.append(
        {
            "domain": "Variance decomposition / transfer stability (1804)",
            "score": s04,
            "status": "PASS" if s04 >= 0.75 else "FAIL",
            "note": d1804.get("verdict", ""),
        }
    )

    global_score = float(sum(float(c["score"]) for c in checks) / len(checks))
    hard_gate = all(c["status"] == "PASS" for c in checks)

    if hard_gate and global_score >= 0.85:
        readiness = "DECOHERENCE_BRANCH_READY"
    elif global_score >= 0.60:
        readiness = "DECOHERENCE_BRANCH_PARTIAL"
    else:
        readiness = "DECOHERENCE_BRANCH_NOT_READY"

    recommendation = (
        "CONTINUE_DECOHERENCE_BRANCH" if hard_gate else "PARK_FOR_PRIMARY_CLAIMS_USE_AS_SECONDARY_DIAGNOSTIC"
    )

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
        "# RAPORT QW-1805: DECOHERENCE BRANCH GATE",
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

    print(f"[QW-1805] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1805] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
