#!/usr/bin/env python3
"""
QW-1819: Sequence branch rigor gate.

Aggregates:
- QW-1815 (in-sample multiscale embedding evidence),
- QW-1816 (regime-aware stabilization attempt),
- QW-1817 (OOS predictive validation),
- QW-1818 (robust OOS attempt).

Outputs a strict readiness verdict and prioritized inconsistency list.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1819_sequence_rigor_gate.json"
OUT_MD = ROOT / "RAPORT_QW1819_SEQUENCE_RIGOR_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def unit_clipped(x: float, lo: float, hi: float) -> float:
    if hi <= lo:
        return 0.0
    return max(0.0, min(1.0, (x - lo) / (hi - lo)))


def main() -> None:
    d1815 = load("report_qw1815_multiscale_sequence_embedding.json")
    d1816 = load("report_qw1816_regime_aware_sequence_model.json")
    d1817 = load("report_qw1817_sequence_oos_validation.json")
    d1818 = load("report_qw1818_robust_sequence_oos.json")

    s15 = d1815["summary"]
    s16 = d1816["summary"]
    s17 = d1817["summary"]
    s18 = d1818["summary"]

    checks: List[Dict[str, object]] = []

    score_insample = 0.0
    score_insample += 0.50 * (1.0 if s15["full_delta_m2e_vs_m2"] > 0.0 else 0.0)
    score_insample += 0.30 * float(s15["prob_m2e_gt_m2"])
    score_insample += 0.20 * float(s15["prob_m2e_gt_flat"])
    checks.append(
        {
            "domain": "In-sample evidence (1815)",
            "score": float(score_insample),
            "status": "PASS" if score_insample >= 0.80 else "FAIL",
            "note": d1815.get("verdict", ""),
        }
    )

    score_oos = 0.0
    score_oos += 0.45 * unit_clipped(float(s17["mean_delta_ll_m2e_vs_m2"]), 0.0, 6.0)
    score_oos += 0.35 * float(s17["prob_delta_ll_m2e_gt_m2"])
    score_oos += 0.20 * float(s17["prob_rmse_gain_positive"])
    checks.append(
        {
            "domain": "OOS predictive support (1817)",
            "score": float(score_oos),
            "status": "PASS" if score_oos >= 0.75 else "FAIL",
            "note": d1817.get("verdict", ""),
        }
    )

    # Stabilization quality: reward lower dispersion and penalize failed stabilization attempts.
    std15 = float(s15["std_delta_m2e_vs_m2"])
    std17 = float(s17["std_delta_ll_m2e_vs_m2"])
    std18 = float(s18["std_delta_ll_m2e_vs_m2"])

    score_stab = 0.0
    score_stab += 0.40 * unit_clipped(0.30 - std15, -0.8, 0.30)  # target <= 0.30
    score_stab += 0.35 * unit_clipped(4.0 - std17, -2.0, 4.0)    # target <= ~4 in OOS LL units
    score_stab += 0.25 * unit_clipped(std17 - std18, -2.0, 2.0)  # improvement in robust attempt
    checks.append(
        {
            "domain": "Stability and dispersion control (1815/1817/1818)",
            "score": float(score_stab),
            "status": "PASS" if score_stab >= 0.65 else "FAIL",
            "note": "dispersion remains high" if score_stab < 0.65 else "controlled",
        }
    )

    # Regime stabilization branch should not degrade vs 1815.
    score_regime = 1.0 if float(s16["full_delta_m2er_vs_m2e"]) > 0 and float(s16["std_reduction_vs_m2e"]) > 0 else 0.0
    checks.append(
        {
            "domain": "Regime-aware stabilization (1816)",
            "score": float(score_regime),
            "status": "PASS" if score_regime >= 0.5 else "FAIL",
            "note": d1816.get("verdict", ""),
        }
    )

    global_score = float(sum(float(c["score"]) for c in checks) / len(checks))
    hard_gate = all(c["status"] == "PASS" for c in checks)

    inconsistencies: List[str] = []
    if float(s15["std_delta_m2e_vs_m2"]) > 0.30:
        inconsistencies.append("In-sample dispersion criterion fails for embedding branch (std delta > 0.30).")
    if float(s17["prob_m2e_gt_flat"]) < 0.95:
        inconsistencies.append("OOS superiority vs flat is not universal (P(M2E>flat) < 0.95).")
    if float(s18["std_reduction_ratio_vs_qw1817"]) <= 0.0:
        inconsistencies.append("Robust preprocessing did not reduce OOS dispersion; variability source remains unresolved.")
    if float(s16["full_delta_m2er_vs_m2e"]) <= 0.0:
        inconsistencies.append("Regime-aware extension degraded evidence vs base embedding model.")

    if hard_gate:
        readiness = "SEQUENCE_BRANCH_READY_FOR_EMPIRICAL_DEPLOYMENT"
        recommendation = "LOCK_PROTOCOL_AND_MOVE_TO_EXTERNAL_VALIDATION"
    elif global_score >= 0.70:
        readiness = "SEQUENCE_BRANCH_PARTIAL_WITH_STABILITY_GAP"
        recommendation = "ADD_HETEROSCEDASTIC_OR_HEAVY_TAIL_LIKELIHOOD_AND_RETEST"
    else:
        readiness = "SEQUENCE_BRANCH_NOT_READY"
        recommendation = "PARK_AND_REDESIGN_MODELING_ASSUMPTIONS"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "global_score": global_score,
        "hard_gate": "PASS" if hard_gate else "FAIL",
        "readiness": readiness,
        "recommendation": recommendation,
        "inconsistencies": inconsistencies,
        "key_metrics": {
            "qw1815_std_delta": std15,
            "qw1817_std_delta_ll": std17,
            "qw1818_std_delta_ll": std18,
            "qw1817_prob_m2e_gt_flat": float(s17["prob_m2e_gt_flat"]),
            "qw1816_full_delta_m2er_vs_m2e": float(s16["full_delta_m2er_vs_m2e"]),
        },
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1819: SEQUENCE RIGOR GATE",
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

    lines.append("")
    lines.append("## Inconsistencies")
    if inconsistencies:
        for it in inconsistencies:
            lines.append(f"- {it}")
    else:
        lines.append("- None identified.")

    lines.extend(["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1819] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1819] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
