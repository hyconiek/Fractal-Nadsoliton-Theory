#!/usr/bin/env python3
"""
QW-1812: Dynamic-proxy branch gate.

Aggregates QW-1808, QW-1809, QW-1811 to decide whether
scalar dynamic proxy models are viable.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1812_dynamic_proxy_branch_gate.json"
OUT_MD = ROOT / "RAPORT_QW1812_DYNAMIC_PROXY_BRANCH_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    d1808 = load("report_qw1808_dynamic_drift_regime_model.json")
    d1809 = load("report_qw1809_dynamic_feature_scan.json")
    d1811 = load("report_qw1811_dynamic_triad_model.json")

    checks: List[Dict[str, object]] = []

    p08 = d1808.get("pass_flags", {})
    s08 = (
        (0.30 if p08.get("full_gain") else 0.0)
        + (0.30 if p08.get("replication_gain") else 0.0)
        + (0.20 if p08.get("dispersion_control") else 0.0)
        + (0.20 if p08.get("dynamic_residual_reduction") else 0.0)
    )
    checks.append(
        {
            "domain": "Drift proxy model (1808)",
            "score": s08,
            "status": "PASS" if s08 >= 0.75 else "FAIL",
            "note": d1808.get("verdict", ""),
        }
    )

    p09 = d1809.get("pass_flags", {})
    s09 = (0.55 if p09.get("gain") else 0.0) + (0.45 if p09.get("dispersion_control") else 0.0)
    checks.append(
        {
            "domain": "Feature scan proxies (1809)",
            "score": s09,
            "status": "PASS" if s09 >= 0.75 else "FAIL",
            "note": d1809.get("verdict", ""),
        }
    )

    p11 = d1811.get("pass_flags", {})
    s11 = (
        (0.35 if p11.get("full_gain") else 0.0)
        + (0.25 if p11.get("replication_gain") else 0.0)
        + (0.20 if p11.get("dispersion_control") else 0.0)
        + (0.20 if p11.get("dynamic_feature_residual_reduction") else 0.0)
    )
    checks.append(
        {
            "domain": "Triad dynamic proxy model (1811)",
            "score": s11,
            "status": "PASS" if s11 >= 0.75 else "FAIL",
            "note": d1811.get("verdict", ""),
        }
    )

    global_score = float(sum(float(c["score"]) for c in checks) / len(checks))
    hard_gate = all(c["status"] == "PASS" for c in checks)

    if hard_gate:
        readiness = "DYNAMIC_PROXY_BRANCH_READY"
        recommendation = "CONTINUE_PROXY_PROGRAM"
    elif global_score >= 0.45:
        readiness = "DYNAMIC_PROXY_BRANCH_PARTIAL"
        recommendation = "USE_PROXIES_ONLY_AS_DIAGNOSTICS"
    else:
        readiness = "DYNAMIC_PROXY_BRANCH_NOT_READY"
        recommendation = "MOVE_TO_SEQUENCE_LEVEL_DYNAMIC_MODELS"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "global_score": global_score,
        "hard_gate": "PASS" if hard_gate else "FAIL",
        "readiness": readiness,
        "recommendation": recommendation,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1812: DYNAMIC PROXY BRANCH GATE",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Global score: {global_score:.3f}",
        f"- Hard gate: **{output['hard_gate']}**",
        f"- Readiness: **{readiness}**",
        f"- Recommendation: **{recommendation}**",
        "",
        "## Checks",
    ]
    for c in checks:
        lines.append(f"- {c['domain']}: {c['status']} | score={c['score']:.3f} | note={c['note']}")
    lines.extend(["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1812] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1812] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
