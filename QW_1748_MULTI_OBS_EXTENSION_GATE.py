#!/usr/bin/env python3
"""
QW-1748: Extension gate after dynamic multi-observable round (1746-1747).
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1748_multi_obs_extension_gate.json"
OUT_MD = ROOT / "RAPORT_QW1748_MULTI_OBS_EXTENSION_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    d1745 = load("report_qw1745_micromodel_rigor_gate.json")
    d1746 = load("report_qw1746_dynamic_observables_derivation.json")
    d1747 = load("report_qw1747_multiobs_joint_inference.json")

    checks: List[Dict[str, object]] = []

    # Baseline carry-over
    base = float(d1745.get("global_score", 0.0))
    checks.append(
        {
            "domain": "Previous rigor baseline (1745)",
            "score": base,
            "status": "PASS" if base >= 0.6 else "FAIL",
            "note": d1745.get("readiness", ""),
        }
    )

    # 1746 scoring
    p46 = d1746.get("pass_flags", {})
    s46 = (0.4 if p46.get("n_good_enough") else 0.0) + (0.3 if p46.get("omega_consistency") else 0.0) + (0.3 if p46.get("spread_control") else 0.0)
    checks.append(
        {
            "domain": "Dynamic observable extraction (1746)",
            "score": s46,
            "status": "PASS" if s46 >= 0.7 else "FAIL",
            "note": d1746.get("verdict", ""),
        }
    )

    # 1747 scoring
    p47 = d1747.get("pass_flags", {})
    s47 = (0.3 if p47.get("conditioning") else 0.0) + (0.25 if p47.get("correlation") else 0.0) + (0.2 if p47.get("bootstrap_width") else 0.0) + (0.25 if p47.get("nonboundary") else 0.0)
    checks.append(
        {
            "domain": "Joint multi-observable identifiability (1747)",
            "score": s47,
            "status": "PASS" if s47 >= 0.7 else "FAIL",
            "note": d1747.get("verdict", ""),
        }
    )

    global_score = float(sum(float(c["score"]) for c in checks) / len(checks))
    hard_gate = all(c["status"] == "PASS" for c in checks)

    if hard_gate and global_score >= 0.8:
        readiness = "MULTI_OBS_EXTENSION_CLOSED"
    elif global_score >= 0.6:
        readiness = "MULTI_OBS_EXTENSION_PARTIAL"
    else:
        readiness = "MULTI_OBS_EXTENSION_OPEN"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "global_score": global_score,
        "hard_gate": "PASS" if hard_gate else "FAIL",
        "readiness": readiness,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1748: MULTI-OBS EXTENSION GATE",
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

    print(f"[QW-1748] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1748] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
