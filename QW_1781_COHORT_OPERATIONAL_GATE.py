#!/usr/bin/env python3
"""
QW-1781: Operational cohort gate for empirical campaign.

Selects operational cohort under pre-registered objective:
1) positive reparam evidence,
2) gain over legacy,
3) Monte-Carlo stability.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1781_cohort_operational_gate.json"
OUT_MD = ROOT / "RAPORT_QW1781_COHORT_OPERATIONAL_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    d1780 = load("report_qw1780_pta_cohort_reparam_bayes.json")
    cohorts = d1780.get("cohort_results", [])

    eligible = [c for c in cohorts if bool(c.get("eligible", False))]
    if len(eligible) == 0:
        output = {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "status": "NO_ELIGIBLE_COHORT",
            "operational_cohort": None,
            "verdict": "COHORT_OPERATIONAL_NOT_READY",
        }
        OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")
        OUT_MD.write_text("# RAPORT QW-1781: COHORT OPERATIONAL GATE\n\n- No eligible cohort.\n", encoding="utf-8")
        print(f"[QW-1781] Saved JSON: {OUT_JSON.name}")
        print(f"[QW-1781] Saved MD:   {OUT_MD.name}")
        return

    # Primary objective for operational selection
    # positive logB first, then delta_logB, then std penalties.
    scored: List[Dict[str, object]] = []
    for c in eligible:
        logb = float(c["logB_reparam_mean"])
        dlog = float(c["delta_logB_mean"])
        s1 = float(c["logB_reparam_std"])
        s2 = float(c["delta_logB_std"])
        pos_bonus = 1.0 if logb > 0.0 else 0.0
        score = 2.0 * pos_bonus + 1.2 * logb + 0.8 * dlog - 1.5 * s1 - 1.0 * s2
        scored.append({"cohort": c, "score": score})

    scored.sort(key=lambda z: float(z["score"]), reverse=True)
    best = scored[0]["cohort"]

    pass_positive = float(best["logB_reparam_mean"]) > 0.0
    pass_gain = float(best["delta_logB_mean"]) >= 0.50
    pass_stability = float(best["logB_reparam_std"]) <= 0.08 and float(best["delta_logB_std"]) <= 0.10
    pass_size = int(best["n_pairs"]) >= 85

    if pass_positive and pass_gain and pass_stability and pass_size:
        verdict = "COHORT_OPERATIONAL_READY"
    elif pass_gain and pass_stability and pass_size:
        verdict = "COHORT_OPERATIONAL_PARTIAL_READY"
    else:
        verdict = "COHORT_OPERATIONAL_NOT_READY"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_eligible": len(eligible),
        "operational_cohort": best,
        "pass_flags": {
            "positive_reparam_logB": bool(pass_positive),
            "gain_over_legacy": bool(pass_gain),
            "mc_stability": bool(pass_stability),
            "enough_pairs": bool(pass_size),
        },
        "verdict": verdict,
        "all_scored": [
            {
                "cohort": s["cohort"]["cohort"],
                "score": float(s["score"]),
                "logB_reparam_mean": float(s["cohort"]["logB_reparam_mean"]),
                "delta_logB_mean": float(s["cohort"]["delta_logB_mean"]),
                "n_pairs": int(s["cohort"]["n_pairs"]),
            }
            for s in scored
        ],
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1781: COHORT OPERATIONAL GATE",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Eligible cohorts: {len(eligible)}",
        f"- Operational cohort: {best['cohort']}",
        f"- logB_reparam_mean: {best['logB_reparam_mean']:.4f}",
        f"- delta_logB_mean: {best['delta_logB_mean']:.4f}",
        f"- n_pairs: {best['n_pairs']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- positive_reparam_logB: {pass_positive}",
        f"- gain_over_legacy: {pass_gain}",
        f"- mc_stability: {pass_stability}",
        f"- enough_pairs: {pass_size}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1781] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1781] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
