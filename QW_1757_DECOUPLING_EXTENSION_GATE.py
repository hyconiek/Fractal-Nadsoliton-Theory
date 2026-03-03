#!/usr/bin/env python3
"""
QW-1757: Gate for decoupling extension (1754-1756).
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1757_decoupling_extension_gate.json"
OUT_MD = ROOT / "RAPORT_QW1757_DECOUPLING_EXTENSION_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    d1754 = load("report_qw1754_orthogonal_extension_gate.json")
    d1755 = load("report_qw1755_beta_null_vs_positive_evidence.json")
    d1756 = load("report_qw1756_omega_beta_decoupling_residual_test.json")

    checks: List[Dict[str, object]] = []

    s54 = float(d1754.get("global_score", 0.0))
    checks.append(
        {
            "domain": "Orthogonal extension baseline (1754)",
            "score": s54,
            "status": "PASS" if s54 >= 0.60 else "FAIL",
            "note": d1754.get("readiness", ""),
        }
    )

    p55 = d1755.get("pass_flags", {})
    s55 = (
        (0.35 if p55.get("evidence_strength") else 0.0)
        + (0.25 if p55.get("positive_beta_probability") else 0.0)
        + (0.20 if p55.get("fit_improvement") else 0.0)
        + (0.20 if p55.get("nonboundary_beta") else 0.0)
    )
    checks.append(
        {
            "domain": "Beta positive-vs-null evidence (1755)",
            "score": s55,
            "status": "PASS" if s55 >= 0.70 else "FAIL",
            "note": d1755.get("verdict", ""),
        }
    )

    p56 = d1756.get("pass_flags", {})
    s56 = (
        (0.15 if p56.get("n_good_enough") else 0.0)
        + (0.20 if p56.get("sensitivity_control") else 0.0)
        + (0.20 if p56.get("orthogonality_proxy_corr") else 0.0)
        + (0.15 if p56.get("spread_control") else 0.0)
        + (0.15 if p56.get("nonboundary_omega") else 0.0)
        + (0.15 if p56.get("fit_quality") else 0.0)
    )
    checks.append(
        {
            "domain": "Omega-beta residual decoupling (1756)",
            "score": s56,
            "status": "PASS" if s56 >= 0.70 else "FAIL",
            "note": d1756.get("verdict", ""),
        }
    )

    global_score = float(sum(float(c["score"]) for c in checks) / len(checks))
    hard_gate = all(c["status"] == "PASS" for c in checks)

    if hard_gate and global_score >= 0.80:
        readiness = "DECOUPLING_EXTENSION_CLOSED"
    elif global_score >= 0.60:
        readiness = "DECOUPLING_EXTENSION_PARTIAL"
    else:
        readiness = "DECOUPLING_EXTENSION_OPEN"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "global_score": global_score,
        "hard_gate": "PASS" if hard_gate else "FAIL",
        "readiness": readiness,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1757: DECOUPLING EXTENSION GATE",
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

    print(f"[QW-1757] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1757] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
