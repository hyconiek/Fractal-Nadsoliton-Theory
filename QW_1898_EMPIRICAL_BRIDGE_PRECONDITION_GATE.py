#!/usr/bin/env python3
"""
QW-1898: Empirical bridge precondition gate.

Aggregates recent corrective branch:
- QW-1896 admissible amplitude-lite single-split pass
- QW-1897 multisplit robustness pass
- complexity admissibility constraints
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict

ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1898_empirical_bridge_precondition_gate.json"
OUT_MD = ROOT / "RAPORT_QW1898_EMPIRICAL_BRIDGE_PRECONDITION_GATE.md"


def read_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def to_score(v: float, lo: float, hi: float) -> float:
    if hi <= lo:
        return 0.0
    x = (float(v) - lo) / (hi - lo)
    return max(0.0, min(1.0, x))


def main() -> None:
    d1896 = read_json("report_qw1896_admissible_amplitude_lite_gate.json")
    d1897 = read_json("report_qw1897_admissible_amplitude_lite_multisplit.json")

    s6 = d1896.get("summaries", {}).get("test", {})
    c6 = d1896.get("complexity", {})
    d6 = d1896.get("delta_vs_1891", {})

    s7 = d1897.get("summary", {})

    score6 = (
        0.35 * to_score(d6.get("rmse_median_gain", 0.0), 0.0, 0.08)
        + 0.25 * to_score(s6.get("canon_median", 0.0), 0.75, 0.98)
        + 0.25 * to_score(s6.get("nonboundary_rate", 0.0), 0.3, 1.0)
        + 0.15 * to_score(c6.get("residual_dof", 0.0), 1.0, 3.0)
    )

    score7 = (
        0.35 * to_score(s7.get("success_rate", 0.0), 0.30, 0.80)
        + 0.25 * to_score(s7.get("rmse_gain_median", 0.0), 0.0, 0.05)
        + 0.20 * to_score(s7.get("canon_gain_median", 0.0), 0.0, 0.10)
        + 0.20 * to_score(s7.get("nonboundary_gain_median", 0.0), 0.0, 0.50)
    )

    global_score = 0.40 * score6 + 0.60 * score7

    hard = bool(
        c6.get("residual_dof", 0) >= 2
        and s7.get("success_rate", 0.0) >= 0.60
        and s7.get("rmse_gain_median", 0.0) > 0.0
        and s7.get("canon_gain_median", 0.0) > 0.0
        and s7.get("nonboundary_gain_median", 0.0) > 0.0
    )

    if hard and global_score >= 0.70:
        readiness = "EMPIRICAL_BRIDGE_PRECONDITION_PASS"
    elif hard:
        readiness = "EMPIRICAL_BRIDGE_PRECONDITION_PARTIAL"
    else:
        readiness = "EMPIRICAL_BRIDGE_PRECONDITION_FAIL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "qw1896": {
                "verdict": d1896.get("verdict", "UNKNOWN"),
                "test_summary": s6,
                "delta_vs_1891": d6,
                "complexity": c6,
                "score": score6,
            },
            "qw1897": {
                "verdict": d1897.get("verdict", "UNKNOWN"),
                "summary": s7,
                "score": score7,
            },
        },
        "global_score": global_score,
        "hard_gate": hard,
        "readiness": readiness,
        "verdict": "EMPIRICAL_BRIDGE_PRECONDITION_GATE_COMPLETE",
        "required_next_step": "QW-1899_EXTERNAL_DETECTOR_PROTOCOL_DESIGN",
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1898: EMPIRICAL BRIDGE PRECONDITION GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        f"- readiness: **{readiness}**",
        f"- hard_gate: **{'PASS' if hard else 'FAIL'}**",
        f"- global_score: {global_score:.3f}",
        "",
        "## Components",
        f"- QW-1896 score: {score6:.3f}",
        f"- QW-1897 score: {score7:.3f}",
        "",
        "## Key Metrics",
        f"- 1897 success_rate: {s7.get('success_rate', 0.0):.3f}",
        f"- 1897 rmse_gain_median: {s7.get('rmse_gain_median', 0.0):.4f}",
        f"- 1897 canon_gain_median: {s7.get('canon_gain_median', 0.0):.4f}",
        f"- 1897 nonboundary_gain_median: {s7.get('nonboundary_gain_median', 0.0):.4f}",
        f"- 1896 residual_dof: {c6.get('residual_dof', 0)}",
        "",
        "## Required Next Step",
        f"- {out['required_next_step']}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1898] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1898] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
