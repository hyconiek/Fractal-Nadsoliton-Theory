#!/usr/bin/env python3
"""
QW-1764: Gate for memory-field iteration (1762 + 1763).
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1764_memory_field_gate.json"
OUT_MD = ROOT / "RAPORT_QW1764_MEMORY_FIELD_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    d1762 = load("report_qw1762_mechanism_augmentation_gate.json")
    d1763 = load("report_qw1763_memory_field_microdynamics_beta_test.json")

    checks: List[Dict[str, object]] = []

    s62 = float(d1762.get("global_score", 0.0))
    checks.append(
        {
            "domain": "Mechanism augmentation baseline (1762)",
            "score": s62,
            "status": "PASS" if s62 >= 0.60 else "FAIL",
            "note": d1762.get("readiness", ""),
        }
    )

    p63 = d1763.get("pass_flags", {})
    s63 = (
        (0.20 if p63.get("n_runs") else 0.0)
        + (0.35 if p63.get("evidence_strength") else 0.0)
        + (0.25 if p63.get("beta_nonboundary") else 0.0)
        + (0.20 if p63.get("fit_improvement") else 0.0)
    )
    checks.append(
        {
            "domain": "Memory-field microdynamics test (1763)",
            "score": s63,
            "status": "PASS" if s63 >= 0.70 else "FAIL",
            "note": d1763.get("verdict", ""),
        }
    )

    global_score = float(sum(float(c["score"]) for c in checks) / len(checks))
    hard_gate = all(c["status"] == "PASS" for c in checks)

    if hard_gate and global_score >= 0.80:
        readiness = "MEMORY_FIELD_ITERATION_CLOSED"
    elif global_score >= 0.60:
        readiness = "MEMORY_FIELD_ITERATION_PARTIAL"
    else:
        readiness = "MEMORY_FIELD_ITERATION_OPEN"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "global_score": global_score,
        "hard_gate": "PASS" if hard_gate else "FAIL",
        "readiness": readiness,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1764: MEMORY FIELD GATE",
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

    print(f"[QW-1764] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1764] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
