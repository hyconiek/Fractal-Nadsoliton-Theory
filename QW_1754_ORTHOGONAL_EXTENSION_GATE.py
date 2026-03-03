#!/usr/bin/env python3
"""
QW-1754: Gate for orthogonal extension strategy (1751-1753).
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1754_orthogonal_extension_gate.json"
OUT_MD = ROOT / "RAPORT_QW1754_ORTHOGONAL_EXTENSION_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    d1751 = load("report_qw1751_beta_separation_gate.json")
    d1752 = load("report_qw1752_omega_orthogonal_observable.json")
    d1753 = load("report_qw1753_orthogonal_triad_joint_inference.json")

    checks: List[Dict[str, object]] = []

    base = float(d1751.get("global_score", 0.0))
    checks.append(
        {
            "domain": "Baseline separation gate (1751)",
            "score": base,
            "status": "PASS" if base >= 0.60 else "FAIL",
            "note": d1751.get("readiness", ""),
        }
    )

    p52 = d1752.get("pass_flags", {})
    s52 = (
        (0.30 if p52.get("fit_quality") else 0.0)
        + (0.25 if p52.get("spread_control") else 0.0)
        + (0.20 if p52.get("beta_orthogonality") else 0.0)
        + (0.25 if p52.get("n_good_enough") else 0.0)
    )
    checks.append(
        {
            "domain": "Omega orthogonal observable (1752)",
            "score": s52,
            "status": "PASS" if s52 >= 0.70 else "FAIL",
            "note": d1752.get("verdict", ""),
        }
    )

    p53 = d1753.get("pass_flags", {})
    s53 = (
        (0.18 if p53.get("beta_observable_quality") else 0.0)
        + (0.18 if p53.get("omega_observable_quality") else 0.0)
        + (0.16 if p53.get("information_content") else 0.0)
        + (0.16 if p53.get("orthogonality") else 0.0)
        + (0.12 if p53.get("phi_stability") else 0.0)
        + (0.20 if p53.get("nonboundary_joint") else 0.0)
    )
    checks.append(
        {
            "domain": "Orthogonal triad joint inference (1753)",
            "score": s53,
            "status": "PASS" if s53 >= 0.70 else "FAIL",
            "note": d1753.get("verdict", ""),
        }
    )

    global_score = float(sum(float(c["score"]) for c in checks) / len(checks))
    hard_gate = all(c["status"] == "PASS" for c in checks)

    if hard_gate and global_score >= 0.80:
        readiness = "ORTHOGONAL_EXTENSION_CLOSED"
    elif global_score >= 0.60:
        readiness = "ORTHOGONAL_EXTENSION_PARTIAL"
    else:
        readiness = "ORTHOGONAL_EXTENSION_OPEN"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "global_score": global_score,
        "hard_gate": "PASS" if hard_gate else "FAIL",
        "readiness": readiness,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1754: ORTHOGONAL EXTENSION GATE",
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

    print(f"[QW-1754] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1754] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
