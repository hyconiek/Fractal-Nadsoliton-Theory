#!/usr/bin/env python3
"""
QW-1919: Stage-B derivational gate after QW-1917/1918.

Integrates:
- Stage-A status from QW-1916,
- no-ansatz triad derivation quality from QW-1917,
- blind external validation of derived triad from QW-1918.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1919_stage_b_derivational_gate.json"
OUT_MD = ROOT / "RAPORT_QW1919_STAGE_B_DERIVATIONAL_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    d1916 = load("report_qw1916_closure_stage_gate.json")
    d1917 = load("report_qw1917_triad_derivation_no_ansatz_profile.json")
    d1918 = load("report_qw1918_triad_blind_external_validation.json")

    stage_a_closed = str(d1916.get("readiness")) == "TOE_STAGE_A_CLOSED_STAGE_B_OPEN"

    v1917 = str(d1917.get("verdict"))
    deriv_strong = v1917 == "TRIAD_NO_ANSATZ_DERIVATION_STRONG"
    deriv_partial = v1917 in {
        "TRIAD_NO_ANSATZ_DERIVATION_STRONG",
        "TRIAD_NO_ANSATZ_DERIVATION_PARTIAL",
    }

    flags1918 = d1918.get("pass_flags", {})
    primary_pass = bool(flags1918.get("primary_all_pass", False))
    stress_soft_pass = bool(flags1918.get("stress_soft_pass", False))

    deriv_strength = 1.0 if deriv_strong else (0.6 if deriv_partial else 0.0)
    score = (
        0.35 * float(stage_a_closed)
        + 0.30 * float(primary_pass)
        + 0.15 * float(stress_soft_pass)
        + 0.20 * float(deriv_strength)
    )

    if stage_a_closed and deriv_strong and primary_pass and stress_soft_pass:
        readiness = "TOE_STAGE_B_CLOSED_PROVISIONAL"
    elif stage_a_closed and primary_pass and deriv_partial:
        readiness = "TOE_STAGE_B_PARTIAL_EXTERNAL_PASS_DERIVATIONAL_PARTIAL"
    elif stage_a_closed and primary_pass:
        readiness = "TOE_STAGE_B_EXTERNAL_PASS_DERIVATIONAL_WEAK"
    else:
        readiness = "TOE_STAGE_B_OPEN"

    if readiness == "TOE_STAGE_B_CLOSED_PROVISIONAL":
        required_next = "PREPARE_FINAL_TEX_INTEGRATION_WITH_FORMAL_LIMITATION_SECTION"
    elif readiness == "TOE_STAGE_B_PARTIAL_EXTERNAL_PASS_DERIVATIONAL_PARTIAL":
        required_next = "RUN_HIGH_POWER_IDENTIFIABILITY_EXPERIMENT_FOR_TRIAD_INTERIOR_STABILITY"
    else:
        required_next = "STRENGTHEN_NO_ANSATZ_IDENTIFIABILITY_AND_REPEAT_BLIND_EXTERNAL_VALIDATION"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "q1916_readiness": d1916.get("readiness"),
            "q1917_verdict": v1917,
            "q1918_verdict": d1918.get("verdict"),
        },
        "flags": {
            "stage_a_closed": bool(stage_a_closed),
            "deriv_strong": bool(deriv_strong),
            "deriv_partial": bool(deriv_partial),
            "primary_blind_external_pass": bool(primary_pass),
            "stress_blind_external_soft_pass": bool(stress_soft_pass),
        },
        "stage_b_score": float(score),
        "readiness": readiness,
        "verdict": "STAGE_B_DERIVATIONAL_GATE_COMPLETE",
        "required_next_step": required_next,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1919: STAGE-B DERIVATIONAL GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        f"- readiness: **{readiness}**",
        f"- stage_b_score: {score:.3f}",
        "",
        "## Flags",
        f"- stage_a_closed: {stage_a_closed}",
        f"- deriv_strong: {deriv_strong}",
        f"- deriv_partial: {deriv_partial}",
        f"- primary_blind_external_pass: {primary_pass}",
        f"- stress_blind_external_soft_pass: {stress_soft_pass}",
        "",
        "## Required Next Step",
        f"- {required_next}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1919] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1919] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1919] readiness={readiness} score={score:.3f}")


if __name__ == "__main__":
    main()
