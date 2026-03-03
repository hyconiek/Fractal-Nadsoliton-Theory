#!/usr/bin/env python3
"""
QW-1851: Joint V2 external confirmatory execution gate.

Combines:
- GW branch readiness (QW-1842 and QW-1839/1841 evidence),
- PTA V2 prereg freeze (QW-1850),
- strict path decision (QW-1849).
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1851_joint_v2_external_execution_gate.json"
OUT_MD = ROOT / "RAPORT_QW1851_JOINT_V2_EXTERNAL_EXECUTION_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def clamp01(x: float) -> float:
    return max(0.0, min(1.0, x))


def main() -> None:
    d1842 = load("report_qw1842_joint_confirmatory_execution_gate.json")
    d1849 = load("report_qw1849_pta_strict_path_selection_gate.json")
    d1850 = load("report_qw1850_pta_v2_prereg_protocol.json")

    pass_gw_operational = d1842.get("hard_gate") == "PASS"
    pass_path_b = d1849.get("decision", {}).get("best_path") == "PATH_B_VERSIONED_CRITERION_REPARAM_WITH_EXTERNAL_CONFIRMATORY"
    pass_pta_v2_frozen = d1850.get("verdict") == "PTA_V2_PREREG_FROZEN_EXTERNAL_CONFIRMATORY_PENDING"

    external_required = bool(d1850.get("protocol", {}).get("status", {}).get("confirmatory_requires_new_external_data", True))

    score = clamp01(
        0.40 * float(pass_gw_operational)
        + 0.25 * float(pass_path_b)
        + 0.35 * float(pass_pta_v2_frozen)
    )

    if pass_gw_operational and pass_path_b and pass_pta_v2_frozen and external_required:
        hard_gate = "PASS"
        readiness = "READY_FOR_EXTERNAL_CONFIRMATORY_V2"
    elif pass_gw_operational and pass_pta_v2_frozen:
        hard_gate = "PARTIAL"
        readiness = "NEAR_READY_PATH_DECISION_MISSING"
    else:
        hard_gate = "FAIL"
        readiness = "NOT_READY_FOR_EXTERNAL_CONFIRMATORY_V2"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "qw1842_hard_gate": d1842.get("hard_gate"),
            "qw1849_best_path": d1849.get("decision", {}).get("best_path"),
            "qw1850_verdict": d1850.get("verdict"),
            "qw1850_protocol_sha256": d1850.get("protocol_sha256"),
        },
        "flags": {
            "gw_operational_ready": bool(pass_gw_operational),
            "pta_path_b_selected": bool(pass_path_b),
            "pta_v2_prereg_frozen": bool(pass_pta_v2_frozen),
            "external_data_required": bool(external_required),
        },
        "strict_score": score,
        "hard_gate": hard_gate,
        "readiness": readiness,
        "verdict": "JOINT_V2_EXTERNAL_EXECUTION_GATE_COMPLETE",
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1851: JOINT V2 EXTERNAL EXECUTION GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Strict score: {score:.3f}",
        f"- Hard gate: **{hard_gate}**",
        f"- Readiness: **{readiness}**",
        "",
        "## Flags",
        f"- gw_operational_ready: {pass_gw_operational}",
        f"- pta_path_b_selected: {pass_path_b}",
        f"- pta_v2_prereg_frozen: {pass_pta_v2_frozen}",
        f"- external_data_required: {external_required}",
        "",
        "## Key Inputs",
        f"- QW-1849 best path: `{d1849.get('decision', {}).get('best_path')}`",
        f"- PTA V2 protocol hash: `{d1850.get('protocol_sha256')}`",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1851] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1851] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
