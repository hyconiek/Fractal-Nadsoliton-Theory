#!/usr/bin/env python3
"""
QW-1876: Kernel closure phase-5 gate.

Aggregates phase-4 with targeted orthogonal forcing and canon-anchored fit.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1876_kernel_closure_phase5_gate.json"
OUT_MD = ROOT / "RAPORT_QW1876_KERNEL_CLOSURE_PHASE5_GATE.md"


def read_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def clip01(x: float) -> float:
    return max(0.0, min(1.0, float(x)))


def main() -> None:
    d1873 = read_json("report_qw1873_kernel_closure_phase4_gate.json")
    d1874 = read_json("report_qw1874_beta_omega_orthogonal_forcing.json")
    d1875 = read_json("report_qw1875_canon_anchored_constrained_fit.json")

    s_prev = clip01(float(d1873.get("global_score", 0.0)))

    s4 = d1874.get("summary", {})
    s_orth = clip01(
        0.25 * float(s4.get("rmse_improved_fraction", 0.0))
        + 0.35 * float(s4.get("corr_reduced_fraction", 0.0))
        + 0.40 * clip01(float(s4.get("orthogonal_canon_median", 0.0)) / (1e-5))
    )

    s5 = d1875.get("summary", {})
    s_anchor = clip01(
        0.30 * float(s5.get("rmse_improved_fraction", 0.0))
        + 0.30 * float(s5.get("canon_improved_fraction", 0.0))
        + 0.40 * clip01(float(s5.get("anchored_canon_median", 0.0)) / (1e-4))
    )
    v1875 = str(d1875.get("verdict", ""))
    if "TRADEOFF_NOT_ACCEPTABLE" in v1875:
        s_anchor *= 0.45

    global_score = 0.40 * s_prev + 0.30 * s_orth + 0.30 * s_anchor

    checks = {
        "phase4_baseline": {
            "score": s_prev,
            "pass": s_prev >= 0.55,
            "note": d1873.get("readiness"),
        },
        "orthogonal_forcing": {
            "score": s_orth,
            "pass": s_orth >= 0.60,
            "note": d1874.get("verdict"),
        },
        "canon_anchored_fit": {
            "score": s_anchor,
            "pass": s_anchor >= 0.60,
            "note": v1875,
        },
    }

    hard_gate = all(v["pass"] for v in checks.values()) and global_score >= 0.66

    if hard_gate:
        readiness = "PHASE5_READY_FOR_EXTERNAL_CONFIRMATORY"
    elif global_score >= 0.48:
        readiness = "PHASE5_PARTIAL_REQUIRES_MODEL_REFORMULATION"
    else:
        readiness = "PHASE5_OPEN_NOT_CLOSED"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "global_score": global_score,
        "hard_gate": hard_gate,
        "readiness": readiness,
        "verdict": "KERNEL_CLOSURE_PHASE5_GATE_COMPLETE",
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1876: KERNEL CLOSURE PHASE-5 GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        f"- Global score: {global_score:.3f}",
        f"- Hard gate: **{'PASS' if hard_gate else 'FAIL'}**",
        f"- Readiness: **{readiness}**",
        "",
        "## Checks",
    ]

    for k, v in checks.items():
        lines.append(f"- {k}: {'PASS' if v['pass'] else 'FAIL'} | score={v['score']:.3f} | note={v['note']}")

    lines += ["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1876] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1876] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
