#!/usr/bin/env python3
"""
QW-1916: Closure stage gate after alpha-bridge integration.

Purpose:
- integrate empirical closure, transfer robustness, and alpha derivational bridge,
- classify current stage of TOE closure in a strict, explicit way.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1916_closure_stage_gate.json"
OUT_MD = ROOT / "RAPORT_QW1916_CLOSURE_STAGE_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    d1902 = load("report_qw1902_empirical_closure_gate.json")
    d1913 = load("report_qw1913_external_pta_multisplit_transfer_stress.json")
    d1915 = load("report_qw1915_alpha_derivational_bridge.json")
    d1890 = load("report_qw1890_toe_closure_decision_gate.json")

    empirical_pass = d1902.get("readiness") == "EMPIRICAL_CLOSURE_PASS"
    transfer_pass = d1913.get("verdict") == "MULTISPLIT_TRANSFER_PASS_ALL_FOLDS"
    alpha_bridge_pass = d1915.get("verdict") == "ALPHA_DERIVATIONAL_BRIDGE_COMPATIBLE"
    derivational_core_closed = bool(d1890.get("hard_gate", False))

    # Stage score: empirical + robustness + bridge + derivational core
    stage_score = (
        0.35 * float(empirical_pass)
        + 0.25 * float(transfer_pass)
        + 0.20 * float(alpha_bridge_pass)
        + 0.20 * float(derivational_core_closed)
    )

    if empirical_pass and transfer_pass and alpha_bridge_pass and derivational_core_closed:
        readiness = "TOE_CLOSURE_READY_STRONG"
    elif empirical_pass and transfer_pass and alpha_bridge_pass and not derivational_core_closed:
        readiness = "TOE_STAGE_A_CLOSED_STAGE_B_OPEN"
    elif empirical_pass and transfer_pass:
        readiness = "TOE_EMPIRICAL_STRONG_BRIDGE_PARTIAL"
    else:
        readiness = "TOE_OPEN"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "q1902_readiness": d1902.get("readiness"),
            "q1913_verdict": d1913.get("verdict"),
            "q1915_verdict": d1915.get("verdict"),
            "q1890_readiness": d1890.get("readiness"),
            "q1890_hard_gate": bool(d1890.get("hard_gate", False)),
        },
        "flags": {
            "empirical_pass": bool(empirical_pass),
            "transfer_pass": bool(transfer_pass),
            "alpha_bridge_pass": bool(alpha_bridge_pass),
            "derivational_core_closed": bool(derivational_core_closed),
        },
        "stage_score": float(stage_score),
        "readiness": readiness,
        "verdict": "CLOSURE_STAGE_GATE_COMPLETE",
        "required_next_step": (
            "DERIVE_BETA_OMEGA_PHI_FROM_MICROMODEL_WITHOUT_ANSATZ_AND_VALIDATE_ON_BLIND_EXTERNAL_DATA"
            if readiness != "TOE_CLOSURE_READY_STRONG"
            else "PREPARE_FINAL_FULL_TOE_PUBLICATION_PACKAGE"
        ),
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1916: CLOSURE STAGE GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        f"- readiness: **{readiness}**",
        f"- stage_score: {stage_score:.3f}",
        "",
        "## Flags",
        f"- empirical_pass: {empirical_pass}",
        f"- transfer_pass: {transfer_pass}",
        f"- alpha_bridge_pass: {alpha_bridge_pass}",
        f"- derivational_core_closed: {derivational_core_closed}",
        "",
        "## Required Next Step",
        f"- {out['required_next_step']}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1916] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1916] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1916] readiness={readiness} stage_score={stage_score:.3f}")


if __name__ == "__main__":
    main()

