#!/usr/bin/env python3
"""
QW-1842: Joint confirmatory execution gate (integrated).

Aggregates:
- QW-1838 global reparam readiness,
- QW-1839 prereg protocol freeze,
- QW-1840 joint power sensitivity,
- QW-1841 GW permutation-null rejection.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1842_joint_confirmatory_execution_gate.json"
OUT_MD = ROOT / "RAPORT_QW1842_JOINT_CONFIRMATORY_EXECUTION_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def clamp01(x: float) -> float:
    return max(0.0, min(1.0, x))


def main() -> None:
    d1838 = load("report_qw1838_global_reparam_readiness_gate.json")
    d1839 = load("report_qw1839_joint_confirmatory_prereg_protocol.json")
    d1840 = load("report_qw1840_joint_confirmatory_power_sensitivity.json")
    d1841 = load("report_qw1841_gw_calibrated_permutation_null.json")

    pass_1838 = d1838.get("hard_gate") == "PASS"
    pass_1839 = d1839.get("verdict") == "JOINT_CONFIRMATORY_PREREG_FROZEN"

    s0 = d1840.get("scenario_results", [])[0]
    p_joint_nom = float(s0.get("pass_probability", {}).get("joint", 0.0)) if s0 else 0.0
    p_joint_ci = s0.get("pass_probability_ci95", {}).get("joint", [0.0, 0.0]) if s0 else [0.0, 0.0]
    pass_1840 = p_joint_nom >= 0.90 and float(p_joint_ci[0]) >= 0.85

    pvals = d1841.get("p_values", {})
    pass_1841 = (
        float(pvals.get("auc_right_tail", 1.0)) <= 0.01
        and float(pvals.get("adv_right_tail", 1.0)) <= 0.01
        and float(pvals.get("gap_left_tail", 1.0)) <= 0.01
    )

    score = clamp01(
        0.30 * float(pass_1838)
        + 0.20 * float(pass_1839)
        + 0.25 * float(pass_1840)
        + 0.25 * float(pass_1841)
    )

    if pass_1838 and pass_1839 and pass_1840 and pass_1841:
        hard_gate = "PASS"
        readiness = "READY_FOR_EXTERNAL_CONFIRMATORY_EXECUTION"
    elif pass_1838 and pass_1839 and (pass_1840 or pass_1841):
        hard_gate = "PARTIAL"
        readiness = "NEAR_READY_CONFIRMATORY_NEEDS_ONE_BLOCK"
    else:
        hard_gate = "FAIL"
        readiness = "NOT_READY_CONFIRMATORY"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "qw1838": {
                "hard_gate": d1838.get("hard_gate"),
                "readiness": d1838.get("readiness"),
            },
            "qw1839": {
                "verdict": d1839.get("verdict"),
                "protocol_sha256": d1839.get("protocol_sha256"),
            },
            "qw1840": {
                "nominal_joint_pass_probability": p_joint_nom,
                "nominal_joint_pass_probability_ci95": p_joint_ci,
            },
            "qw1841": {
                "p_values": pvals,
                "verdict": d1841.get("verdict"),
            },
        },
        "pass_flags": {
            "global_reparam_ready": bool(pass_1838),
            "prereg_frozen": bool(pass_1839),
            "joint_power_sufficient": bool(pass_1840),
            "gw_null_rejected": bool(pass_1841),
        },
        "global_score": score,
        "hard_gate": hard_gate,
        "readiness": readiness,
        "verdict": "JOINT_CONFIRMATORY_EXECUTION_GATE_COMPLETE",
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1842: JOINT CONFIRMATORY EXECUTION GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Global score: {score:.3f}",
        f"- Hard gate: **{hard_gate}**",
        f"- Readiness: **{readiness}**",
        "",
        "## Pass Flags",
        f"- global_reparam_ready (1838): {pass_1838}",
        f"- prereg_frozen (1839): {pass_1839}",
        f"- joint_power_sufficient (1840): {pass_1840}",
        f"- gw_null_rejected (1841): {pass_1841}",
        "",
        "## Key Inputs",
        f"- QW-1839 protocol hash: `{d1839.get('protocol_sha256')}`",
        f"- QW-1840 nominal joint pass prob: {p_joint_nom:.3f} (CI95: [{float(p_joint_ci[0]):.3f}, {float(p_joint_ci[1]):.3f}])",
        (
            "- QW-1841 p-values (AUC/ADV/GAP): "
            f"{float(pvals.get('auc_right_tail', 1.0)):.6f}, "
            f"{float(pvals.get('adv_right_tail', 1.0)):.6f}, "
            f"{float(pvals.get('gap_left_tail', 1.0)):.6f}"
        ),
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1842] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1842] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
