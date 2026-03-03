#!/usr/bin/env python3
"""
QW-1873: Kernel closure phase-4 gate.

Aggregates:
- Phase-3B status (QW-1868),
- Primary node evidence corpus quality (QW-1871),
- Structural node dynamic bridge test (QW-1872).
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1873_kernel_closure_phase4_gate.json"
OUT_MD = ROOT / "RAPORT_QW1873_KERNEL_CLOSURE_PHASE4_GATE.md"


def read_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def clip01(x: float) -> float:
    return max(0.0, min(1.0, float(x)))


def main() -> None:
    d1868 = read_json("report_qw1868_kernel_closure_phase3b_gate.json")
    d1871 = read_json("report_qw1871_primary_node_evidence_corpus.json")
    d1872 = read_json("report_qw1872_structural_node_dynamic_micromodel.json")

    s_prev = clip01(float(d1868.get("global_score", 0.0)))

    sum71 = d1871.get("summary", {})
    n_primary = float(sum71.get("n_primary_files_with_node_hit", 0.0))
    p_win = float(d1871.get("winner_posterior", 0.0))
    s_node_data = clip01(0.55 * clip01(n_primary / 15.0) + 0.45 * p_win)

    sum72 = d1872.get("summary", {})
    canon_med = float(sum72.get("node_canon_median", 0.0))
    canon_gain = float(sum72.get("median_canon_gain", 0.0))
    rmse_frac = float(sum72.get("rmse_improved_fraction", 0.0))
    c_imp = float(sum72.get("canon_improved_fraction", 0.0))

    s_bridge = clip01(
        0.35 * clip01(canon_med / 0.25)
        + 0.25 * clip01((canon_gain + 0.02) / 0.10)
        + 0.20 * rmse_frac
        + 0.20 * c_imp
    )

    global_score = 0.45 * s_prev + 0.20 * s_node_data + 0.35 * s_bridge

    checks = {
        "phase3b_status": {
            "score": s_prev,
            "pass": s_prev >= 0.50,
            "note": d1868.get("readiness"),
        },
        "primary_node_data_quality": {
            "score": s_node_data,
            "pass": s_node_data >= 0.45,
            "note": d1871.get("verdict"),
        },
        "structural_bridge": {
            "score": s_bridge,
            "pass": s_bridge >= 0.55,
            "note": d1872.get("verdict"),
        },
    }

    hard_gate = all(v["pass"] for v in checks.values()) and global_score >= 0.62

    if hard_gate:
        readiness = "PHASE4_READY_FOR_EXTERNAL_KERNEL_CONFIRMATORY"
    elif global_score >= 0.45:
        readiness = "PHASE4_PARTIAL_NODE_DATA_COLLECTION_REQUIRED"
    else:
        readiness = "PHASE4_OPEN_STRUCTURAL_REDESIGN_REQUIRED"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "global_score": global_score,
        "hard_gate": hard_gate,
        "readiness": readiness,
        "verdict": "KERNEL_CLOSURE_PHASE4_GATE_COMPLETE",
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1873: KERNEL CLOSURE PHASE-4 GATE",
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

    print(f"[QW-1873] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1873] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
