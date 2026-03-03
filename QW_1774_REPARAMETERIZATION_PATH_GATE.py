#!/usr/bin/env python3
"""
QW-1774: Reparameterization path gate.

Aggregates:
- QW-1768 (leakage-controlled non-envelope mechanism support)
- QW-1772 (kernel bridge integration baseline)
- QW-1773 (omega-suppressed projection to legacy scale)
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1774_reparameterization_path_gate.json"
OUT_MD = ROOT / "RAPORT_QW1774_REPARAMETERIZATION_PATH_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    d1768 = load("report_qw1768_leakage_controlled_paired_intervention.json")
    d1772 = load("report_qw1772_kernel_bridge_integration_gate.json")
    d1773 = load("report_qw1773_omega_suppressed_legacy_projection.json")

    checks: List[Dict[str, object]] = []

    # 1768: mechanism branch validity
    p68 = d1768.get("pass_flags", {})
    s68 = (
        (0.30 if p68.get("discrimination") else 0.0)
        + (0.30 if p68.get("continuous_regression") else 0.0)
        + (0.20 if p68.get("calibration_nonboundary") else 0.0)
        + (0.20 if p68.get("cv_stability") else 0.0)
    )
    checks.append(
        {
            "domain": "Leakage-controlled mechanism branch (1768)",
            "score": s68,
            "status": "PASS" if s68 >= 0.70 else "FAIL",
            "note": d1768.get("verdict", ""),
        }
    )

    # 1772: prior bridge baseline
    s72 = float(d1772.get("global_score", 0.0))
    checks.append(
        {
            "domain": "Pre-reparameterization bridge baseline (1772)",
            "score": s72,
            "status": "PASS" if s72 >= 0.60 else "FAIL",
            "note": d1772.get("readiness", ""),
        }
    )

    # 1773: new projection path quality
    p73 = d1773.get("pass_flags", {})
    s73 = (
        (0.35 if p73.get("distance_reduction") else 0.0)
        + (0.25 if p73.get("ci_overlap") else 0.0)
        + (0.20 if p73.get("parsimony") else 0.0)
        + (0.20 if p73.get("projection_shape") else 0.0)
    )
    checks.append(
        {
            "domain": "Omega-suppressed projection quality (1773)",
            "score": s73,
            "status": "PASS" if s73 >= 0.70 else "FAIL",
            "note": d1773.get("verdict", ""),
        }
    )

    global_score = float(sum(float(c["score"]) for c in checks) / len(checks))
    hard_gate = all(c["status"] == "PASS" for c in checks)

    pass_reparam_core = bool(p73.get("distance_reduction")) and bool(p73.get("ci_overlap")) and bool(p73.get("parsimony"))
    pass_mech_core = bool(p68.get("discrimination")) and bool(p68.get("continuous_regression"))

    if hard_gate and global_score >= 0.80:
        readiness = "REPARAMETERIZATION_PATH_CLOSED"
    elif global_score >= 0.70 and pass_reparam_core and pass_mech_core:
        readiness = "REPARAMETERIZATION_PATH_STRONG_PARTIAL"
    elif global_score >= 0.55:
        readiness = "REPARAMETERIZATION_PATH_PARTIAL"
    else:
        readiness = "REPARAMETERIZATION_PATH_OPEN"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "global_score": global_score,
        "hard_gate": "PASS" if hard_gate else "FAIL",
        "readiness": readiness,
        "core_flags": {
            "mechanism_branch_core": pass_mech_core,
            "reparameterization_core": pass_reparam_core,
        },
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1774: REPARAMETERIZATION PATH GATE",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Global score: {global_score:.3f}",
        f"- Hard gate: **{output['hard_gate']}**",
        f"- Readiness: **{readiness}**",
        "",
        "## Core Flags",
        f"- mechanism_branch_core: {pass_mech_core}",
        f"- reparameterization_core: {pass_reparam_core}",
        "",
        "## Checks",
    ]
    for c in checks:
        lines.append(f"- {c['domain']}: {c['status']} | score={c['score']:.3f} | note={c['note']}")
    lines.extend(["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1774] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1774] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
