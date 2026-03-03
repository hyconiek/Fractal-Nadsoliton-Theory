#!/usr/bin/env python3
"""
QW-1890: Global closure decision gate after signed-coupling branch.

Aggregates phase-7 node-state status (QW-1886) with signed-coupling
multisplit holdout tuning (QW-1889).
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict

ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1890_toe_closure_decision_gate.json"
OUT_MD = ROOT / "RAPORT_QW1890_TOE_CLOSURE_DECISION_GATE.md"


def read_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def to_score(v: float, lo: float, hi: float) -> float:
    if hi <= lo:
        return 0.0
    x = (float(v) - lo) / (hi - lo)
    return max(0.0, min(1.0, x))


def main() -> None:
    d1886 = read_json("report_qw1886_node_state_phase7_gate.json")
    d1889 = read_json("report_qw1889_signed_coupling_multisplit_tuning.json")

    phase7_score = float(d1886.get("global_score", 0.0))
    phase7_hard = bool(d1886.get("hard_gate", False))

    hold = d1889.get("holdout_summary", {})
    succ = float(hold.get("success_rate", 0.0))
    rmse_gain = float(hold.get("rmse_gain_median", 0.0))
    canon_gain = float(hold.get("canon_gain_median", -1.0))
    canon_q25 = float(hold.get("canon_gain_q25", -1.0))
    nb_gain = float(hold.get("nonboundary_gain_median", 0.0))

    signed_score = (
        0.35 * to_score(succ, 0.20, 0.70)
        + 0.25 * to_score(rmse_gain, 0.0, 0.03)
        + 0.20 * to_score(nb_gain, 0.0, 0.50)
        + 0.20 * to_score(canon_gain, -0.10, 0.05)
    )

    signed_hard = bool(succ >= 0.50 and canon_gain >= -0.05 and canon_q25 >= -0.20)

    global_score = 0.60 * phase7_score + 0.40 * signed_score
    hard_gate = bool(phase7_hard and signed_hard)

    if hard_gate and global_score >= 0.70:
        readiness = "TOE_CLOSURE_READY_FOR_EMPIRICAL_BRIDGE"
    elif signed_hard and not phase7_hard:
        readiness = "TOE_PARTIAL_PROGRESS_REQUIRES_STRONGER_THEORETICAL_CONSTRAINTS"
    else:
        readiness = "TOE_NOT_CLOSED_REQUIRES_DERIVATIONAL_REFORMULATION"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "qw1886": {
                "global_score": phase7_score,
                "hard_gate": phase7_hard,
                "readiness": d1886.get("readiness", "UNKNOWN"),
            },
            "qw1889": {
                "holdout_summary": hold,
                "signed_score": signed_score,
                "signed_hard_gate": signed_hard,
                "verdict": d1889.get("verdict", "UNKNOWN"),
            },
        },
        "global_score": global_score,
        "hard_gate": hard_gate,
        "readiness": readiness,
        "verdict": "TOE_CLOSURE_DECISION_GATE_COMPLETE",
        "required_next_step": "QW-1891_DERIVATIONAL_CONSTRAINTS_FROM_NADSOLITON",
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1890: TOE CLOSURE DECISION GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        f"- readiness: **{readiness}**",
        f"- hard_gate: **{'PASS' if hard_gate else 'FAIL'}**",
        f"- global_score: {global_score:.3f}",
        "",
        "## Components",
        f"- QW-1886 phase7: score={phase7_score:.3f}, hard={phase7_hard}",
        f"- QW-1889 signed-holdout: score={signed_score:.3f}, hard={signed_hard}",
        "",
        "## QW-1889 Key Holdout Metrics",
        f"- success_rate: {succ:.3f}",
        f"- rmse_gain_median: {rmse_gain:.4f}",
        f"- canon_gain_median: {canon_gain:.4f}",
        f"- canon_gain_q25: {canon_q25:.4f}",
        f"- nonboundary_gain_median: {nb_gain:.4f}",
        "",
        "## Required Next Step",
        f"- {out['required_next_step']}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1890] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1890] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
