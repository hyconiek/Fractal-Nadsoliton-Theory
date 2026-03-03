#!/usr/bin/env python3
"""
QW-1878: Mechanism comparison (legacy vs reformulated core).

Compares QW-1874 (legacy orthogonal forcing) vs QW-1877 (node-state reform)
under preregistered metrics.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1878_reformulation_mechanism_comparison.json"
OUT_MD = ROOT / "RAPORT_QW1878_REFORMULATION_MECHANISM_COMPARISON.md"


def read_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def clip01(x: float) -> float:
    return max(0.0, min(1.0, float(x)))


def main() -> None:
    d1874 = read_json("report_qw1874_beta_omega_orthogonal_forcing.json")
    d1877 = read_json("report_qw1877_node_state_structural_reformulation.json")

    s4 = d1874.get("summary", {})
    s7 = d1877.get("summary", {})

    delta = {
        "rmse_median_change": float(s7.get("reform_rmse_median", 0.0)) - float(s4.get("orthogonal_rmse_median", 0.0)),
        "canon_median_change": float(s7.get("reform_canon_median", 0.0)) - float(s4.get("orthogonal_canon_median", 0.0)),
        "rmse_improved_fraction_change": float(s7.get("rmse_improved_fraction", 0.0)) - float(s4.get("rmse_improved_fraction", 0.0)),
        "canon_improved_fraction_change": float(s7.get("canon_improved_fraction", 0.0)) - float(s4.get("canon_improved_fraction", 0.0)),
    }

    score_rmse = clip01(0.5 + 2.0 * (-delta["rmse_median_change"]))
    score_canon = clip01(0.5 + 2e5 * delta["canon_median_change"])
    score_frac = clip01(0.5 + 0.5 * delta["rmse_improved_fraction_change"] + 0.5 * delta["canon_improved_fraction_change"])

    comparison_score = 0.35 * score_rmse + 0.40 * score_canon + 0.25 * score_frac

    if comparison_score >= 0.65:
        verdict = "REFORMULATION_MECHANISM_SUPERIOR"
    elif comparison_score >= 0.50:
        verdict = "REFORMULATION_MECHANISM_MIXED"
    else:
        verdict = "REFORMULATION_MECHANISM_NOT_SUPERIOR"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "legacy_source": "QW-1874",
        "reform_source": "QW-1877",
        "delta": delta,
        "subscores": {
            "rmse": score_rmse,
            "canon": score_canon,
            "fractions": score_frac,
        },
        "comparison_score": comparison_score,
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1878: REFORMULATION MECHANISM COMPARISON",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- comparison score: {comparison_score:.3f}",
        "",
        "## Delta (QW-1877 minus QW-1874)",
        f"- rmse median change: {delta['rmse_median_change']:.4e}",
        f"- canon median change: {delta['canon_median_change']:.4e}",
        f"- rmse improved fraction change: {delta['rmse_improved_fraction_change']:.4e}",
        f"- canon improved fraction change: {delta['canon_improved_fraction_change']:.4e}",
        "",
        "## Subscores",
        f"- rmse: {score_rmse:.3f}",
        f"- canon: {score_canon:.3f}",
        f"- fractions: {score_frac:.3f}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1878] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1878] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
