#!/usr/bin/env python3
"""
QW-1793: Operational model lock-in gate.

Integrates results from 1786-1792 to lock one default empirical protocol.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1793_model_lockin_gate.json"
OUT_MD = ROOT / "RAPORT_QW1793_MODEL_LOCKIN_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    d1786 = load("report_qw1786_campaign_robustness_gate.json")
    d1788 = load("report_qw1788_q_prior_width_sweep.json")
    d1789 = load("report_qw1789_cohort_criteria_sweep.json")
    d1790 = load("report_qw1790_residual_angular_augmentation.json")
    d1791 = load("report_qw1791_robust_likelihood_reparam.json")
    d1792 = load("report_qw1792_heteroscedastic_reparam.json")

    checks: List[Dict[str, object]] = []

    # 1) Campaign robustness must hold.
    s86 = float(d1786.get("global_score", 0.0))
    pass86 = d1786.get("hard_gate") == "PASS" and s86 >= 0.90
    checks.append(
        {
            "domain": "Campaign robustness (1786)",
            "score": s86,
            "status": "PASS" if pass86 else "FAIL",
            "note": d1786.get("readiness", ""),
        }
    )

    # 2) Hyperparameter lock (q-width and cohort criteria) must be strong.
    pass88 = d1788.get("verdict") == "Q_PRIOR_WIDTH_SELECTION_SUPPORTED" and d1788.get("recommendation_strength") == "STRONG"
    pass89 = d1789.get("verdict") == "COHORT_CRITERIA_SELECTION_SUPPORTED" and d1789.get("recommendation_strength") == "STRONG"
    s_hyper = (0.5 if pass88 else 0.0) + (0.5 if pass89 else 0.0)
    checks.append(
        {
            "domain": "Hyperparameter lock (1788+1789)",
            "score": s_hyper,
            "status": "PASS" if s_hyper >= 1.0 else "FAIL",
            "note": f"q_width={d1788.get('recommended_q_width')}, criteria={d1789.get('recommended_criteria', {}).get('name')}",
        }
    )

    # 3) Angular augmentation must be rejected (negative control).
    s90 = d1790.get("summary", {})
    delta90 = float(s90.get("full_delta_logB_m3_minus_m2", 0.0))
    pass90 = d1790.get("verdict") == "ANGULAR_AUGMENTATION_WEAK" and delta90 < 0.0
    checks.append(
        {
            "domain": "Negative-control extension rejected (1790)",
            "score": 1.0 if pass90 else 0.0,
            "status": "PASS" if pass90 else "FAIL",
            "note": f"delta_M3-M2={delta90:.4f}",
        }
    )

    # 4) Alternative likelihoods should not beat baseline in absolute logB.
    s91 = d1791.get("summary", {})
    s92 = d1792.get("summary", {})
    delta_t = float(s91.get("full_delta_t_minus_gauss", 0.0))
    delta_h = float(s92.get("full_gain_hetero_minus_homo_reparam_vs_flat", 0.0))
    pass91 = delta_t <= 0.0
    pass92 = delta_h <= 0.0
    s_like = (0.5 if pass91 else 0.0) + (0.5 if pass92 else 0.0)
    checks.append(
        {
            "domain": "Likelihood alternatives vs baseline (1791+1792)",
            "score": s_like,
            "status": "PASS" if s_like >= 1.0 else "FAIL",
            "note": f"delta_t-gauss={delta_t:.4f}, delta_hetero-homo={delta_h:.4f}",
        }
    )

    global_score = float(sum(float(c["score"]) for c in checks) / len(checks))
    hard_gate = all(c["status"] == "PASS" for c in checks)

    if hard_gate and global_score >= 0.95:
        readiness = "MODEL_LOCKIN_CONFIRMED"
    elif global_score >= 0.75:
        readiness = "MODEL_LOCKIN_PARTIAL"
    else:
        readiness = "MODEL_LOCKIN_NOT_CONFIRMED"

    operational_protocol = {
        "signal_model": "M2_reparam_HDq_plus_constant",
        "noise_model": "homoscedastic_gaussian",
        "fraction": 0.95,
        "q_width": float(d1788.get("recommended_q_width", 0.20)),
        "cohort": d1789.get("recommended_criteria", {}),
        "frozen_after_qw": 1793,
    }

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "global_score": global_score,
        "hard_gate": "PASS" if hard_gate else "FAIL",
        "readiness": readiness,
        "operational_protocol": operational_protocol,
        "model_competitors_rejected": {
            "angular_P2_augmentation_qw1790": {
                "delta_M3_minus_M2": delta90,
                "verdict": d1790.get("verdict"),
            },
            "student_t_likelihood_qw1791": {
                "delta_t_minus_gauss": delta_t,
                "verdict": d1791.get("verdict"),
            },
            "heteroscedastic_likelihood_qw1792": {
                "delta_hetero_minus_homo": delta_h,
                "verdict": d1792.get("verdict"),
            },
        },
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1793: MODEL LOCK-IN GATE",
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
            "## Locked Operational Protocol",
            f"- signal_model: {operational_protocol['signal_model']}",
            f"- noise_model: {operational_protocol['noise_model']}",
            f"- fraction: {operational_protocol['fraction']}",
            f"- q_width: {operational_protocol['q_width']}",
            f"- cohort: {operational_protocol['cohort'].get('name')}",
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1793] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1793] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
