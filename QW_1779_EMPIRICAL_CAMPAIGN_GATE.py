#!/usr/bin/env python3
"""
QW-1779: Empirical campaign gate after first reparam PTA run.

Aggregates:
- QW-1776 readiness
- QW-1777 precheck
- QW-1778 first empirical Bayes run
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1779_empirical_campaign_gate.json"
OUT_MD = ROOT / "RAPORT_QW1779_EMPIRICAL_CAMPAIGN_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    d1776 = load("report_qw1776_reparam_empirical_readiness_gate.json")
    d1777 = load("report_qw1777_empirical_pta_reparam_precheck.json")
    d1778 = load("report_qw1778_pta_bayes_reparam_reanalysis.json")

    checks: List[Dict[str, object]] = []

    s76 = float(d1776.get("global_score", 0.0))
    checks.append(
        {
            "domain": "Readiness gate (1776)",
            "score": s76,
            "status": "PASS" if s76 >= 0.80 and d1776.get("hard_gate", "FAIL") == "PASS" else "FAIL",
            "note": d1776.get("readiness", ""),
        }
    )

    p77 = d1777.get("pass_flags", {})
    s77 = (
        (0.20 if p77.get("reparam_ready") else 0.0)
        + (0.20 if p77.get("pipeline_sensitivity") else 0.0)
        + (0.20 if p77.get("real_structure_signal") else 0.0)
        + (0.20 if p77.get("legacy_path_not_robust") else 0.0)
        + (0.20 if p77.get("bayes_headroom_for_reparam") else 0.0)
    )
    checks.append(
        {
            "domain": "Empirical precheck (1777)",
            "score": s77,
            "status": "PASS" if s77 >= 0.80 else "FAIL",
            "note": d1777.get("verdict", ""),
        }
    )

    p78 = d1778.get("pass_flags", {})
    s78 = (
        (0.45 if p78.get("reparam_gain_over_legacy") else 0.0)
        + (0.25 if p78.get("reparam_nonnegative_evidence") else 0.0)
        + (0.30 if p78.get("enough_valid_pairs") else 0.0)
    )
    checks.append(
        {
            "domain": "First empirical Bayes run (1778)",
            "score": s78,
            "status": "PASS" if s78 >= 0.75 else "FAIL",
            "note": d1778.get("verdict", ""),
        }
    )

    global_score = float(sum(float(c["score"]) for c in checks) / len(checks))
    hard_gate = all(c["status"] == "PASS" for c in checks)

    if hard_gate and global_score >= 0.85:
        readiness = "EMPIRICAL_CAMPAIGN_STRONGLY_ON_TRACK"
    elif global_score >= 0.70:
        readiness = "EMPIRICAL_CAMPAIGN_PARTIAL_ON_TRACK"
    else:
        readiness = "EMPIRICAL_CAMPAIGN_NOT_ON_TRACK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "global_score": global_score,
        "hard_gate": "PASS" if hard_gate else "FAIL",
        "readiness": readiness,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1779: EMPIRICAL CAMPAIGN GATE",
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

    print(f"[QW-1779] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1779] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
