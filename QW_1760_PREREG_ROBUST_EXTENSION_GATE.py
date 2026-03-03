#!/usr/bin/env python3
"""
QW-1760: Gate for prereg + robust extension (1757-1759).
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1760_prereg_robust_extension_gate.json"
OUT_MD = ROOT / "RAPORT_QW1760_PREREG_ROBUST_EXTENSION_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def score_1758(pass_flags: Dict[str, bool]) -> float:
    return (
        (0.12 if pass_flags.get("n_good_enough") else 0.0)
        + (0.16 if pass_flags.get("beta_evidence") else 0.0)
        + (0.12 if pass_flags.get("beta_nonboundary") else 0.0)
        + (0.12 if pass_flags.get("omega_nonboundary") else 0.0)
        + (0.10 if pass_flags.get("beta_spread") else 0.0)
        + (0.10 if pass_flags.get("omega_spread") else 0.0)
        + (0.10 if pass_flags.get("beta_omega_coupling_low") else 0.0)
        + (0.08 if pass_flags.get("beta_perturb_stability") else 0.0)
        + (0.10 if pass_flags.get("fit_quality") else 0.0)
    )


def score_1759(pass_flags: Dict[str, bool]) -> float:
    return (
        (0.12 if pass_flags.get("n_good_enough") else 0.0)
        + (0.16 if pass_flags.get("beta_evidence") else 0.0)
        + (0.12 if pass_flags.get("beta_nonboundary") else 0.0)
        + (0.12 if pass_flags.get("omega_nonboundary") else 0.0)
        + (0.10 if pass_flags.get("beta_spread") else 0.0)
        + (0.10 if pass_flags.get("omega_spread") else 0.0)
        + (0.10 if pass_flags.get("beta_omega_coupling_low") else 0.0)
        + (0.08 if pass_flags.get("beta_perturb_stability") else 0.0)
        + (0.10 if pass_flags.get("phase_quality") else 0.0)
    )


def main() -> None:
    d1757 = load("report_qw1757_decoupling_extension_gate.json")
    d1758 = load("report_qw1758_preregistered_identifiability_protocol.json")
    d1759 = load("report_qw1759_robust_phase_prereg_protocol.json")

    checks: List[Dict[str, object]] = []

    s57 = float(d1757.get("global_score", 0.0))
    checks.append(
        {
            "domain": "Decoupling baseline gate (1757)",
            "score": s57,
            "status": "PASS" if s57 >= 0.60 else "FAIL",
            "note": d1757.get("readiness", ""),
        }
    )

    p58 = d1758.get("pass_flags", {})
    s58 = float(score_1758(p58))
    checks.append(
        {
            "domain": "Strict prereg protocol (1758)",
            "score": s58,
            "status": "PASS" if s58 >= 0.70 else "FAIL",
            "note": d1758.get("verdict", ""),
        }
    )

    p59 = d1759.get("pass_flags", {})
    s59 = float(score_1759(p59))
    checks.append(
        {
            "domain": "Robust-phase prereg protocol (1759)",
            "score": s59,
            "status": "PASS" if s59 >= 0.70 else "FAIL",
            "note": d1759.get("verdict", ""),
        }
    )

    global_score = float(sum(float(c["score"]) for c in checks) / len(checks))
    hard_gate = all(c["status"] == "PASS" for c in checks)

    if hard_gate and global_score >= 0.80:
        readiness = "PREREG_ROBUST_EXTENSION_CLOSED"
    elif global_score >= 0.60:
        readiness = "PREREG_ROBUST_EXTENSION_PARTIAL"
    else:
        readiness = "PREREG_ROBUST_EXTENSION_OPEN"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "global_score": global_score,
        "hard_gate": "PASS" if hard_gate else "FAIL",
        "readiness": readiness,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1760: PREREG ROBUST EXTENSION GATE",
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

    print(f"[QW-1760] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1760] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
