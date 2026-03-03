#!/usr/bin/env python3
"""
QW-1786: Empirical campaign robustness recovery gate.

Integrates:
- campaign launch gate (1779),
- operational cohort readiness (1781),
- low-coverage stress replications (1782, 1784),
- high-coverage stratified replication (1785).
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1786_campaign_robustness_gate.json"
OUT_MD = ROOT / "RAPORT_QW1786_CAMPAIGN_ROBUSTNESS_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def clip01(x: float) -> float:
    return max(0.0, min(1.0, float(x)))


def stress_protocol_score(d: Dict) -> float:
    s = d.get("summary", {})
    full_reparam = float(s.get("full_logB_reparam", 0.0))
    full_delta = float(s.get("full_delta_logB", 0.0))
    p_reparam = float(s.get("prob_logB_reparam_positive", 0.0))
    p_delta = float(s.get("prob_delta_logB_positive", 0.0))

    score = 0.0
    score += 0.30 if full_reparam > 0.0 else 0.0
    score += 0.30 if full_delta > 0.0 else 0.0
    score += 0.20 * clip01(p_reparam / 0.80)
    score += 0.20 * clip01(p_delta / 0.95)
    return float(score)


def main() -> None:
    d1779 = load("report_qw1779_empirical_campaign_gate.json")
    d1781 = load("report_qw1781_cohort_operational_gate.json")
    d1782 = load("report_qw1782_operational_cohort_replication.json")
    d1784 = load("report_qw1784_stratified_replication.json")
    d1785 = load("report_qw1785_high_coverage_replication.json")

    checks: List[Dict[str, object]] = []

    s79 = float(d1779.get("global_score", 0.0))
    pass79 = d1779.get("hard_gate") == "PASS" and s79 >= 0.80
    checks.append(
        {
            "domain": "Campaign launch status (1779)",
            "score": s79,
            "status": "PASS" if pass79 else "FAIL",
            "note": d1779.get("readiness", ""),
        }
    )

    p81 = d1781.get("pass_flags", {})
    s81 = (
        (0.25 if p81.get("positive_reparam_logB") else 0.0)
        + (0.30 if p81.get("gain_over_legacy") else 0.0)
        + (0.25 if p81.get("mc_stability") else 0.0)
        + (0.20 if p81.get("enough_pairs") else 0.0)
    )
    pass81 = s81 >= 0.75
    checks.append(
        {
            "domain": "Operational cohort readiness (1781)",
            "score": s81,
            "status": "PASS" if pass81 else "FAIL",
            "note": d1781.get("verdict", ""),
        }
    )

    s82 = stress_protocol_score(d1782)
    s84 = stress_protocol_score(d1784)
    s_stress = 0.5 * (s82 + s84)
    pass_stress = s_stress >= 0.85
    checks.append(
        {
            "domain": "Low-coverage stress robustness (1782+1784)",
            "score": s_stress,
            "status": "PASS" if pass_stress else "FAIL",
            "note": f"score_1782={s82:.3f}, score_1784={s84:.3f}",
        }
    )

    p85 = d1785.get("pass_flags", {})
    s85 = (
        (0.20 if p85.get("full_cohort_positive") else 0.0)
        + (0.30 if p85.get("positive_reparam_probability") else 0.0)
        + (0.30 if p85.get("gain_probability") else 0.0)
        + (0.20 if p85.get("dispersion_control") else 0.0)
    )
    pass85 = s85 >= 0.90 and d1785.get("verdict") == "HIGH_COVERAGE_REPLICATION_SUPPORTED"
    checks.append(
        {
            "domain": "High-coverage replication stability (1785)",
            "score": s85,
            "status": "PASS" if pass85 else "FAIL",
            "note": d1785.get("verdict", ""),
        }
    )

    global_score = float(sum(float(c["score"]) for c in checks) / len(checks))
    hard_gate = pass79 and pass81 and pass_stress and pass85

    if hard_gate and global_score >= 0.90:
        readiness = "EMPIRICAL_CAMPAIGN_ROBUSTNESS_CONFIRMED"
    elif global_score >= 0.75:
        readiness = "EMPIRICAL_CAMPAIGN_ROBUSTNESS_PARTIAL"
    else:
        readiness = "EMPIRICAL_CAMPAIGN_ROBUSTNESS_NOT_CONFIRMED"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "global_score": global_score,
        "hard_gate": "PASS" if hard_gate else "FAIL",
        "readiness": readiness,
        "stress_protocol_details": {
            "score_qw1782": s82,
            "score_qw1784": s84,
            "mean_stress_score": s_stress,
        },
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1786: CAMPAIGN ROBUSTNESS GATE",
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
    lines.extend(
        [
            "",
            "## Stress Details",
            f"- score_1782: {s82:.3f}",
            f"- score_1784: {s84:.3f}",
            f"- mean_stress_score: {s_stress:.3f}",
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1786] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1786] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
