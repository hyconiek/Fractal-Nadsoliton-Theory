#!/usr/bin/env python3
"""
QW-1864: Phase-3 gate for kernel closure readiness.

Aggregates:
- QW-1861 canonical reconstruction strength,
- QW-1862 micromodel compatibility,
- QW-1863 identifiability design quality.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1864_kernel_closure_phase3_gate.json"
OUT_MD = ROOT / "RAPORT_QW1864_KERNEL_CLOSURE_PHASE3_GATE.md"


def read_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def clip01(x: float) -> float:
    return max(0.0, min(1.0, float(x)))


def main() -> None:
    d1861 = read_json("report_qw1861_canonical_kernel_reconstruction_700_1600.json")
    d1862 = read_json("report_qw1862_micromodel_canonical_compatibility_audit.json")
    d1863 = read_json("report_qw1863_identifiability_optimal_observables_design.json")

    s1861 = clip01(
        float(d1861.get("consistency", {}).get("min_consensus", 0.0))
        * float(d1861.get("coverage", {}).get("evidence_coverage_fraction", 0.0))
    )

    s1862 = clip01(float(d1862.get("summary", {}).get("compatibility_index", 0.0)))

    best_inter_score = float(d1863.get("best_with_signed_intervention", {}).get("score", 0.0))
    cond_gain = float(d1863.get("intervention_gain", {}).get("cond_q90_reduction_factor", 0.0))
    s1863 = clip01(0.65 * clip01(best_inter_score / 0.25) + 0.35 * clip01(cond_gain / 2.0))

    global_score = 0.40 * s1861 + 0.35 * s1862 + 0.25 * s1863

    checks = {
        "canonical_reconstruction_strength": {
            "score": s1861,
            "pass": s1861 >= 0.50,
            "note": d1861.get("verdict"),
        },
        "micromodel_compatibility": {
            "score": s1862,
            "pass": s1862 >= 0.50,
            "note": d1862.get("verdict"),
        },
        "identifiability_design": {
            "score": s1863,
            "pass": s1863 >= 0.55,
            "note": d1863.get("verdict"),
        },
    }

    hard_gate = all(v["pass"] for v in checks.values()) and global_score >= 0.65

    if hard_gate:
        readiness = "PHASE3_READY_FOR_EXTERNAL_CONFIRMATORY_EXECUTION"
    elif global_score >= 0.45:
        readiness = "PHASE3_PARTIAL_READY_SYNTHETIC_RECOVERY_REQUIRED"
    else:
        readiness = "PHASE3_OPEN_KERNEL_NOT_CLOSED"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "global_score": global_score,
        "hard_gate": hard_gate,
        "readiness": readiness,
        "recommended_next_steps": d1863.get("recommended_next_studies", []),
        "verdict": "KERNEL_CLOSURE_PHASE3_GATE_COMPLETE",
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1864: KERNEL CLOSURE PHASE-3 GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        f"- Global score: {global_score:.3f}",
        f"- Hard gate: **{'PASS' if hard_gate else 'FAIL'}**",
        f"- Readiness: **{readiness}**",
        "",
        "## Checks",
    ]

    for name, row in checks.items():
        lines.append(
            f"- {name}: {'PASS' if row['pass'] else 'FAIL'} | score={row['score']:.3f} | note={row['note']}"
        )

    lines += ["", "## Recommended Next Steps"]
    if not out["recommended_next_steps"]:
        lines.append("- none")
    else:
        for s in out["recommended_next_steps"]:
            lines.append(f"- {s.get('id')}: {s.get('title')}")

    lines += ["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1864] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1864] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
