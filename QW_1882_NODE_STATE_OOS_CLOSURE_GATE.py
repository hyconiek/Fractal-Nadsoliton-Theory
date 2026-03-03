#!/usr/bin/env python3
"""
QW-1882: OOS closure gate for node-state branch.

Aggregates strict OOS execution (QW-1880) and identifiability stress (QW-1881).
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1882_node_state_oos_closure_gate.json"
OUT_MD = ROOT / "RAPORT_QW1882_NODE_STATE_OOS_CLOSURE_GATE.md"


def read_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def clip01(x: float) -> float:
    return max(0.0, min(1.0, float(x)))


def main() -> None:
    d1880 = read_json("report_qw1880_node_state_strict_oos.json")
    d1881 = read_json("report_qw1881_node_state_identifiability_stress.json")

    s0 = d1880.get("summaries", {})
    test = s0.get("test", {})
    test_base = s0.get("test_baseline", {})

    oos_score = clip01(
        0.35 * clip01(float(test.get("canon_median", 0.0)) / 0.25)
        + 0.20 * clip01(float(test.get("nonboundary_rate", 0.0)))
        + 0.20 * clip01((float(test_base.get("rmse_median", 1.0)) - float(test.get("rmse_median", 1.0))) / 0.03 + 0.5)
        + 0.25 * clip01((float(test.get("canon_median", 0.0)) - float(test_base.get("canon_median", 0.0))) / 0.20 + 0.5)
    )
    v1880 = str(d1880.get("verdict", ""))
    if "FAIL" in v1880:
        oos_score *= 0.60

    st = d1881.get("global_metrics", {})
    stress_score = clip01(
        0.40 * float(st.get("stability_index", 0.0))
        + 0.20 * clip01(float(st.get("nonboundary_rate", 0.0)))
        + 0.20 * clip01(1.0 - abs(float(st.get("corr_omega_beta", 1.0))) / 0.7)
        + 0.20 * clip01(1.0 - float(st.get("beta_iqr", 1.0)) / 0.08)
    )
    v1881 = str(d1881.get("verdict", ""))
    if "FAIL" in v1881:
        stress_score *= 0.70

    global_score = 0.55 * oos_score + 0.45 * stress_score

    checks = {
        "strict_oos": {
            "score": oos_score,
            "pass": oos_score >= 0.60,
            "note": v1880,
        },
        "identifiability_stress": {
            "score": stress_score,
            "pass": stress_score >= 0.60,
            "note": v1881,
        },
    }

    hard_gate = all(v["pass"] for v in checks.values()) and global_score >= 0.65

    if hard_gate:
        readiness = "NODE_STATE_OOS_CLOSURE_PASS"
    elif global_score >= 0.50:
        readiness = "NODE_STATE_OOS_PARTIAL_REQUIRES_TARGETED_REFINEMENT"
    else:
        readiness = "NODE_STATE_OOS_FAIL_REFORMULATION_NOT_STABLE"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "global_score": global_score,
        "hard_gate": hard_gate,
        "readiness": readiness,
        "verdict": "NODE_STATE_OOS_CLOSURE_GATE_COMPLETE",
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1882: NODE-STATE OOS CLOSURE GATE",
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

    print(f"[QW-1882] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1882] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
