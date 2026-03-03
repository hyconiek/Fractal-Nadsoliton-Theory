#!/usr/bin/env python3
"""
QW-1783: Campaign stability gate after operational replication.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1783_campaign_stability_gate.json"
OUT_MD = ROOT / "RAPORT_QW1783_CAMPAIGN_STABILITY_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    d1779 = load("report_qw1779_empirical_campaign_gate.json")
    d1781 = load("report_qw1781_cohort_operational_gate.json")
    d1782 = load("report_qw1782_operational_cohort_replication.json")

    checks: List[Dict[str, object]] = []

    s79 = float(d1779.get("global_score", 0.0))
    checks.append(
        {
            "domain": "Campaign launch status (1779)",
            "score": s79,
            "status": "PASS" if s79 >= 0.80 and d1779.get("hard_gate") == "PASS" else "FAIL",
            "note": d1779.get("readiness", ""),
        }
    )

    p81 = d1781.get("pass_flags", {})
    s81 = (
        (0.25 if p81.get("positive_reparam_logB") else 0.0)
        + (0.30 if p81.get("gain_over_legacy") else 0.0)
        + (0.25 if p81.get("mc_stability") else 0.0)
        + (0.20 if p81.get("enough_pairs") else 0.0)
    )
    checks.append(
        {
            "domain": "Operational cohort selection (1781)",
            "score": s81,
            "status": "PASS" if s81 >= 0.75 else "FAIL",
            "note": d1781.get("verdict", ""),
        }
    )

    p82 = d1782.get("pass_flags", {})
    s82 = (
        (0.35 if p82.get("positive_reparam_probability") else 0.0)
        + (0.35 if p82.get("gain_probability") else 0.0)
        + (0.30 if p82.get("dispersion_control") else 0.0)
    )
    checks.append(
        {
            "domain": "Operational replication stability (1782)",
            "score": s82,
            "status": "PASS" if s82 >= 0.75 else "FAIL",
            "note": d1782.get("verdict", ""),
        }
    )

    global_score = float(sum(float(c["score"]) for c in checks) / len(checks))
    hard_gate = all(c["status"] == "PASS" for c in checks)

    if hard_gate and global_score >= 0.82:
        readiness = "CAMPAIGN_STABILITY_CONFIRMED"
    elif global_score >= 0.65:
        readiness = "CAMPAIGN_STABILITY_PARTIAL"
    else:
        readiness = "CAMPAIGN_STABILITY_NOT_CONFIRMED"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "global_score": global_score,
        "hard_gate": "PASS" if hard_gate else "FAIL",
        "readiness": readiness,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1783: CAMPAIGN STABILITY GATE",
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

    print(f"[QW-1783] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1783] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
