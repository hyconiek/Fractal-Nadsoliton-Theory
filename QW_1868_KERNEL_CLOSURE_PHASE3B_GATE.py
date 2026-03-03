#!/usr/bin/env python3
"""
QW-1868: Kernel closure phase-3B gate after adaptive benchmarks.

Aggregates phase-3 gate with executed synthetic recovery and adaptive
beta-augmentation benchmarks.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1868_kernel_closure_phase3b_gate.json"
OUT_MD = ROOT / "RAPORT_QW1868_KERNEL_CLOSURE_PHASE3B_GATE.md"


def read_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def clip01(x: float) -> float:
    return max(0.0, min(1.0, float(x)))


def main() -> None:
    d1864 = read_json("report_qw1864_kernel_closure_phase3_gate.json")
    d1865 = read_json("report_qw1865_synthetic_recovery_benchmark.json")
    d1866 = read_json("report_qw1866_paired_signed_intervention_benchmark.json")
    d1867 = read_json("report_qw1867_beta_augmented_observables_benchmark.json")

    checks1864 = d1864.get("checks", {})
    s_theory = float(checks1864.get("canonical_reconstruction_strength", {}).get("score", 0.0))
    s_comp = float(checks1864.get("micromodel_compatibility", {}).get("score", 0.0))
    s_design = float(checks1864.get("identifiability_design", {}).get("score", 0.0))

    m1865 = d1865.get("metrics", {})
    tol5 = m1865.get("tolerance_hit_rate", {})
    s_synth = clip01(
        0.20 * float(tol5.get("omega_0p08", 0.0))
        + 0.20 * float(tol5.get("phi_0p20", 0.0))
        + 0.40 * float(tol5.get("beta_0p015", 0.0))
        + 0.20 * float(m1865.get("nonboundary_rate", 0.0))
    )

    # Penalize if paired signed intervention is clearly negative.
    imp1866 = d1866.get("improvement", {})
    paired_penalty = 1.0
    if float(imp1866.get("rmse_beta_factor", 1.0)) < 1.0:
        paired_penalty *= 0.80
    if float(imp1866.get("tol_beta_gain", 0.0)) < 0.0:
        paired_penalty *= 0.85

    best1867 = d1867.get("best_arm", {})
    s_adapt = clip01(
        0.60 * clip01((float(best1867.get("beta_rmse_factor", 1.0)) - 1.0) / 0.30)
        + 0.40 * clip01((float(best1867.get("beta_tol_gain", 0.0)) + 0.02) / 0.12)
    )
    s_adapt *= paired_penalty

    global_score = (
        0.25 * s_theory
        + 0.25 * s_comp
        + 0.15 * s_design
        + 0.20 * s_synth
        + 0.15 * s_adapt
    )

    checks = {
        "theory_canonical_strength": {
            "score": s_theory,
            "pass": s_theory >= 0.50,
            "note": checks1864.get("canonical_reconstruction_strength", {}).get("note"),
        },
        "micromodel_canonical_compatibility": {
            "score": s_comp,
            "pass": s_comp >= 0.50,
            "note": checks1864.get("micromodel_compatibility", {}).get("note"),
        },
        "design_identifiability": {
            "score": s_design,
            "pass": s_design >= 0.55,
            "note": checks1864.get("identifiability_design", {}).get("note"),
        },
        "synthetic_recovery_execution": {
            "score": s_synth,
            "pass": s_synth >= 0.70,
            "note": d1865.get("verdict"),
        },
        "adaptive_beta_augmentation": {
            "score": s_adapt,
            "pass": s_adapt >= 0.55,
            "note": d1867.get("verdict"),
        },
    }

    hard_gate = all(v["pass"] for v in checks.values()) and global_score >= 0.65

    if hard_gate:
        readiness = "PHASE3B_READY_FOR_REAL_DATA_CONFIRMATORY_CAMPAIGN"
    elif global_score >= 0.45:
        readiness = "PHASE3B_PARTIAL_MICROMODEL_REDESIGN_REQUIRED"
    else:
        readiness = "PHASE3B_OPEN_NOT_CLOSED"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "paired_penalty": paired_penalty,
        "global_score": global_score,
        "hard_gate": hard_gate,
        "readiness": readiness,
        "verdict": "KERNEL_CLOSURE_PHASE3B_GATE_COMPLETE",
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1868: KERNEL CLOSURE PHASE-3B GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        f"- Global score: {global_score:.3f}",
        f"- Hard gate: **{'PASS' if hard_gate else 'FAIL'}**",
        f"- Readiness: **{readiness}**",
        f"- Paired penalty factor: {paired_penalty:.3f}",
        "",
        "## Checks",
    ]

    for name, row in checks.items():
        lines.append(
            f"- {name}: {'PASS' if row['pass'] else 'FAIL'} | score={row['score']:.3f} | note={row['note']}"
        )

    lines += [
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1868] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1868] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
