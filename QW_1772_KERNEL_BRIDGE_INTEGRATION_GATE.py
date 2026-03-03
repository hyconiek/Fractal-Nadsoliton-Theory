#!/usr/bin/env python3
"""
QW-1772: Kernel-bridge integration gate (1769 + 1770 + 1771).
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1772_kernel_bridge_integration_gate.json"
OUT_MD = ROOT / "RAPORT_QW1772_KERNEL_BRIDGE_INTEGRATION_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    d1769 = load("report_qw1769_nonenvelope_rigor_gate.json")
    d1770 = load("report_qw1770_kernel_bridge_from_nonenvelope.json")
    d1771 = load("report_qw1771_bridge_legacy_consistency.json")

    checks: List[Dict[str, object]] = []

    s69 = float(d1769.get("global_score", 0.0))
    checks.append(
        {
            "domain": "Non-envelope branch baseline (1769)",
            "score": s69,
            "status": "PASS" if s69 >= 0.65 else "FAIL",
            "note": d1769.get("readiness", ""),
        }
    )

    p70 = d1770.get("pass_flags", {})
    s70 = (
        (0.24 if p70.get("beta_bridge_regression") else 0.0)
        + (0.24 if p70.get("omega_regression") else 0.0)
        + (0.18 if p70.get("phi_regression") else 0.0)
        + (0.17 if p70.get("nonboundary_predictions") else 0.0)
        + (0.17 if p70.get("cv_stability") else 0.0)
    )
    checks.append(
        {
            "domain": "Kernel bridge regression quality (1770)",
            "score": s70,
            "status": "PASS" if s70 >= 0.70 else "FAIL",
            "note": d1770.get("verdict", ""),
        }
    )

    p71 = d1771.get("pass_flags", {})
    s71 = (
        (0.35 if p71.get("direct_consistency") else 0.0)
        + (0.30 if p71.get("transformed_consistency") else 0.0)
        + (0.20 if p71.get("simple_map") else 0.0)
        + (0.15 if p71.get("monotonicity") else 0.0)
    )
    checks.append(
        {
            "domain": "Bridge-legacy scale consistency (1771)",
            "score": s71,
            "status": "PASS" if s71 >= 0.70 else "FAIL",
            "note": d1771.get("verdict", ""),
        }
    )

    global_score = float(sum(float(c["score"]) for c in checks) / len(checks))
    hard_gate = all(c["status"] == "PASS" for c in checks)

    if hard_gate and global_score >= 0.80:
        readiness = "KERNEL_BRIDGE_INTEGRATION_CLOSED"
    elif global_score >= 0.60:
        readiness = "KERNEL_BRIDGE_INTEGRATION_PARTIAL"
    else:
        readiness = "KERNEL_BRIDGE_INTEGRATION_OPEN"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "global_score": global_score,
        "hard_gate": "PASS" if hard_gate else "FAIL",
        "readiness": readiness,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1772: KERNEL BRIDGE INTEGRATION GATE",
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

    print(f"[QW-1772] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1772] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
