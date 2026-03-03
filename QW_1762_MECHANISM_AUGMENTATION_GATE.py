#!/usr/bin/env python3
"""
QW-1762: Gate for mechanism augmentation attempt (1760 + 1761).
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1762_mechanism_augmentation_gate.json"
OUT_MD = ROOT / "RAPORT_QW1762_MECHANISM_AUGMENTATION_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    d1760 = load("report_qw1760_prereg_robust_extension_gate.json")
    d1761 = load("report_qw1761_envelope_hybrid_model_selection.json")

    checks: List[Dict[str, object]] = []

    s60 = float(d1760.get("global_score", 0.0))
    checks.append(
        {
            "domain": "Prereg+robust baseline (1760)",
            "score": s60,
            "status": "PASS" if s60 >= 0.60 else "FAIL",
            "note": d1760.get("readiness", ""),
        }
    )

    p61 = d1761.get("pass_flags", {})
    s61 = (
        (0.25 if p61.get("hybrid_preferred_rate") else 0.0)
        + (0.20 if p61.get("hybrid_gain_vs_m1") else 0.0)
        + (0.20 if p61.get("hybrid_gain_vs_m0") else 0.0)
        + (0.20 if p61.get("m2_beta_nonboundary") else 0.0)
        + (0.15 if p61.get("m2_lambda_positive") else 0.0)
    )
    checks.append(
        {
            "domain": "Envelope mechanism augmentation (1761)",
            "score": s61,
            "status": "PASS" if s61 >= 0.70 else "FAIL",
            "note": d1761.get("verdict", ""),
        }
    )

    global_score = float(sum(float(c["score"]) for c in checks) / len(checks))
    hard_gate = all(c["status"] == "PASS" for c in checks)

    if hard_gate and global_score >= 0.80:
        readiness = "MECHANISM_AUGMENTATION_CLOSED"
    elif global_score >= 0.60:
        readiness = "MECHANISM_AUGMENTATION_PARTIAL"
    else:
        readiness = "MECHANISM_AUGMENTATION_OPEN"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "global_score": global_score,
        "hard_gate": "PASS" if hard_gate else "FAIL",
        "readiness": readiness,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1762: MECHANISM AUGMENTATION GATE",
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

    print(f"[QW-1762] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1762] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
