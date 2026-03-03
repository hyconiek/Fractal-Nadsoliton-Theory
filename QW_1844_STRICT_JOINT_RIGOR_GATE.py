#!/usr/bin/env python3
"""
QW-1844: Strict joint rigor gate.

Integrates QW-1842 execution readiness with QW-1843 inferential rigor for PTA.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1844_strict_joint_rigor_gate.json"
OUT_MD = ROOT / "RAPORT_QW1844_STRICT_JOINT_RIGOR_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    d1842 = load("report_qw1842_joint_confirmatory_execution_gate.json")
    d1843 = load("report_qw1843_pta_threshold_inference_rigor.json")

    pass_exec = d1842.get("hard_gate") == "PASS"
    pta_flags = d1843.get("pass_flags", {})

    pass_pta_mean = bool(pta_flags.get("mean_threshold_with_ci", False))
    pass_pta_std = bool(pta_flags.get("std_threshold_with_ci", False))
    pass_pta_prob_inf = bool(pta_flags.get("prob_threshold_inferential", False))

    n_now = int(d1843.get("n_replications", 0))
    n_min = int(d1843.get("probability_inference", {}).get("n_min_all_positive_for_alpha0p05_vs_threshold", 0))
    n_additional_if_all_positive = max(0, n_min - n_now)

    score = (
        0.40 * float(pass_exec)
        + 0.25 * float(pass_pta_mean)
        + 0.15 * float(pass_pta_std)
        + 0.20 * float(pass_pta_prob_inf)
    )

    if pass_exec and pass_pta_mean and pass_pta_std and pass_pta_prob_inf:
        hard_gate = "PASS"
        readiness = "STRICT_CONFIRMATORY_READY"
    elif pass_exec and pass_pta_mean and pass_pta_std and (not pass_pta_prob_inf):
        hard_gate = "PARTIAL"
        readiness = "PTA_PROBABILITY_INFERENCE_UNDERPOWERED"
    else:
        hard_gate = "FAIL"
        readiness = "STRICT_NOT_READY"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "qw1842_hard_gate": d1842.get("hard_gate"),
            "qw1842_readiness": d1842.get("readiness"),
            "qw1843_verdict": d1843.get("verdict"),
        },
        "pass_flags": {
            "execution_gate_pass": bool(pass_exec),
            "pta_mean_threshold_with_ci": bool(pass_pta_mean),
            "pta_std_threshold_with_ci": bool(pass_pta_std),
            "pta_prob_threshold_inferential": bool(pass_pta_prob_inf),
        },
        "power_gap": {
            "n_replications_current": n_now,
            "n_replications_min_if_all_positive": n_min,
            "additional_replications_if_all_positive": n_additional_if_all_positive,
        },
        "global_score": float(score),
        "hard_gate": hard_gate,
        "readiness": readiness,
        "verdict": "STRICT_JOINT_RIGOR_GATE_COMPLETE",
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1844: STRICT JOINT RIGOR GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Global score: {score:.3f}",
        f"- Hard gate: **{hard_gate}**",
        f"- Readiness: **{readiness}**",
        "",
        "## Pass Flags",
        f"- execution_gate_pass (1842): {pass_exec}",
        f"- pta_mean_threshold_with_ci: {pass_pta_mean}",
        f"- pta_std_threshold_with_ci: {pass_pta_std}",
        f"- pta_prob_threshold_inferential: {pass_pta_prob_inf}",
        "",
        "## Power Gap (PTA prob threshold)",
        f"- n current: {n_now}",
        f"- n min (all-positive) for alpha=0.05: {n_min}",
        f"- additional replications needed (all-positive case): {n_additional_if_all_positive}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1844] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1844] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
