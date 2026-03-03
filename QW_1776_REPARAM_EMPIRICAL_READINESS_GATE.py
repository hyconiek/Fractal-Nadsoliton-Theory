#!/usr/bin/env python3
"""
QW-1776: Reparameterization empirical-readiness gate.

Aggregates:
- QW-1768 (mechanism branch rigor)
- QW-1774 (reparameterization path gate)
- QW-1775 (OOD validation)
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1776_reparam_empirical_readiness_gate.json"
OUT_MD = ROOT / "RAPORT_QW1776_REPARAM_EMPIRICAL_READINESS_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    d1768 = load("report_qw1768_leakage_controlled_paired_intervention.json")
    d1774 = load("report_qw1774_reparameterization_path_gate.json")
    d1775 = load("report_qw1775_reparam_ood_validation.json")

    checks: List[Dict[str, object]] = []

    p68 = d1768.get("pass_flags", {})
    s68 = (
        (0.30 if p68.get("discrimination") else 0.0)
        + (0.30 if p68.get("continuous_regression") else 0.0)
        + (0.20 if p68.get("calibration_nonboundary") else 0.0)
        + (0.20 if p68.get("cv_stability") else 0.0)
    )
    checks.append(
        {
            "domain": "Mechanism branch rigor (1768)",
            "score": s68,
            "status": "PASS" if s68 >= 0.70 else "FAIL",
            "note": d1768.get("verdict", ""),
        }
    )

    s74 = float(d1774.get("global_score", 0.0))
    checks.append(
        {
            "domain": "Reparameterization path strength (1774)",
            "score": s74,
            "status": "PASS" if s74 >= 0.70 else "FAIL",
            "note": d1774.get("readiness", ""),
        }
    )

    p75 = d1775.get("pass_flags", {})
    s75 = (
        (0.45 if p75.get("core_alignment") else 0.0)
        + (0.25 if p75.get("projection_shape") else 0.0)
        + (0.30 if p75.get("param_sensitivity") else 0.0)
    )
    checks.append(
        {
            "domain": "OOD validation of reparameterization (1775)",
            "score": s75,
            "status": "PASS" if s75 >= 0.70 else "FAIL",
            "note": d1775.get("verdict", ""),
        }
    )

    global_score = float(sum(float(c["score"]) for c in checks) / len(checks))
    hard_gate = all(c["status"] == "PASS" for c in checks)

    if hard_gate and global_score >= 0.82:
        readiness = "REPARAM_EMPIRICAL_READY"
    elif global_score >= 0.70:
        readiness = "REPARAM_PATH_STRONG_PARTIAL_READY"
    else:
        readiness = "REPARAM_PATH_NOT_READY"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "global_score": global_score,
        "hard_gate": "PASS" if hard_gate else "FAIL",
        "readiness": readiness,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1776: REPARAM EMPIRICAL READINESS GATE",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Global score: {global_score:.3f}",
        f"- Hard gate: **{output['hard_gate']}**",
        f"- Readiness: **{readiness}**",
        "",
        "## Checks",
    ]
    for c in checks:
        lines.append(f"- {c['domain']}: {c['status']} | score={c['score']:.3f} | note={c['note']}")
    lines.extend(["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1776] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1776] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
