#!/usr/bin/env python3
"""
QW-1751: Gate for beta-separation strategy (1749-1750).
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1751_beta_separation_gate.json"
OUT_MD = ROOT / "RAPORT_QW1751_BETA_SEPARATION_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    d1748 = load("report_qw1748_multi_obs_extension_gate.json")
    d1749 = load("report_qw1749_beta_orthogonal_observable.json")
    d1750 = load("report_qw1750_separated_beta_then_omega_phi.json")

    checks: List[Dict[str, object]] = []

    # Carry-over baseline
    base = float(d1748.get("global_score", 0.0))
    checks.append(
        {
            "domain": "Pre-separation baseline (1748)",
            "score": base,
            "status": "PASS" if base >= 0.6 else "FAIL",
            "note": d1748.get("readiness", ""),
        }
    )

    # 1749
    p49 = d1749.get("pass_flags", {})
    s49 = (0.3 if p49.get("fit_quality") else 0.0) + (0.25 if p49.get("spread_control") else 0.0) + (0.25 if p49.get("phase_orthogonality") else 0.0) + (0.20 if p49.get("nonboundary_beta") else 0.0)
    checks.append(
        {
            "domain": "Orthogonal beta observable (1749)",
            "score": s49,
            "status": "PASS" if s49 >= 0.7 else "FAIL",
            "note": d1749.get("verdict", ""),
        }
    )

    # 1750
    p50 = d1750.get("pass_flags", {})
    s50 = (0.20 if p50.get("beta_quality") else 0.0) + (0.15 if p50.get("omega_width") else 0.0) + (0.15 if p50.get("phi_width") else 0.0) + (0.20 if p50.get("conditioning_or_corr") else 0.0) + (0.15 if p50.get("nonboundary_omega") else 0.0) + (0.15 if p50.get("beta_sweep_nonboundary") else 0.0)
    checks.append(
        {
            "domain": "Separated omega/phi inference (1750)",
            "score": s50,
            "status": "PASS" if s50 >= 0.7 else "FAIL",
            "note": d1750.get("verdict", ""),
        }
    )

    global_score = float(sum(float(c["score"]) for c in checks) / len(checks))
    hard_gate = all(c["status"] == "PASS" for c in checks)

    if hard_gate and global_score >= 0.8:
        readiness = "BETA_SEPARATION_STRATEGY_CLOSED"
    elif global_score >= 0.6:
        readiness = "BETA_SEPARATION_STRATEGY_PARTIAL"
    else:
        readiness = "BETA_SEPARATION_STRATEGY_OPEN"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "global_score": global_score,
        "hard_gate": "PASS" if hard_gate else "FAIL",
        "readiness": readiness,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1751: BETA SEPARATION GATE",
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

    print(f"[QW-1751] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1751] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
