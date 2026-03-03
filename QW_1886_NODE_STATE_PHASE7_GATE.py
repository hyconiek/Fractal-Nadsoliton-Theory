#!/usr/bin/env python3
"""
QW-1886: Phase-7 gate for node-state closure readiness.

Aggregates:
- QW-1882 strict OOS gate
- QW-1884 Pareto rebalancing
- QW-1885 multisplit robustness

Purpose: decide if current node-state branch is scientifically closure-ready.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict

ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1886_node_state_phase7_gate.json"
OUT_MD = ROOT / "RAPORT_QW1886_NODE_STATE_PHASE7_GATE.md"


def read_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def to_score(v: float, lo: float = 0.0, hi: float = 1.0) -> float:
    if hi <= lo:
        return 0.0
    x = (float(v) - lo) / (hi - lo)
    return max(0.0, min(1.0, x))


def main() -> None:
    d1882 = read_json("report_qw1882_node_state_oos_closure_gate.json")
    d1884 = read_json("report_qw1884_node_state_pareto_oos_rebalancing.json")
    d1885 = read_json("report_qw1885_node_state_multisplit_tradeoff_robustness.json")

    # 1882 strict OOS component.
    s1882 = d1882.get("global_score", 0.0)
    strict_pass = bool(d1882.get("hard_gate", False))

    # 1884 tradeoff component.
    s4 = d1884.get("summaries", {}).get("test", {})
    rmse_ratio = float(d1884.get("delta_vs_1880", {}).get("rmse_ratio", 9.9))
    canon = float(s4.get("canon_median", 0.0))
    nb = float(s4.get("nonboundary_rate", 0.0))

    pareto_score = (
        0.40 * to_score(nb, 0.0, 0.7)
        + 0.30 * to_score(canon, 0.2, 0.95)
        + 0.30 * to_score(1.60 - rmse_ratio, 0.0, 0.60)
    )
    pareto_soft_pass = bool(nb >= 0.50 and canon >= 0.80 and rmse_ratio <= 1.60)

    # 1885 robustness component.
    s5 = d1885.get("summary", {})
    succ = float(s5.get("success_rate", 0.0))
    rmse_pen = float(s5.get("rmse_penalty_median", 9.9))
    nb_gain = float(s5.get("nonboundary_gain_median", 0.0))

    robust_score = (
        0.45 * to_score(succ, 0.10, 0.70)
        + 0.30 * to_score(nb_gain, 0.10, 0.50)
        + 0.25 * to_score(0.10 - rmse_pen, -0.02, 0.10)
    )
    robust_pass = bool(succ >= 0.60 and nb_gain >= 0.30 and rmse_pen <= 0.07)

    global_score = 0.35 * float(s1882) + 0.25 * pareto_score + 0.40 * robust_score

    hard_gate = bool(strict_pass and robust_pass)

    if hard_gate and global_score >= 0.70:
        readiness = "PHASE7_NODE_STATE_CLOSURE_READY"
    elif pareto_soft_pass and not robust_pass:
        readiness = "PHASE7_PARTIAL_SINGLE_SPLIT_ONLY_REQUIRES_MODEL_REFORMULATION"
    else:
        readiness = "PHASE7_NOT_CLOSED_REQUIRES_MICRODYNAMICS_REBUILD"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "qw1882": {
                "hard_gate": strict_pass,
                "score": float(s1882),
                "verdict": d1882.get("verdict", "UNKNOWN"),
            },
            "qw1884": {
                "pareto_score": pareto_score,
                "nonboundary_rate": nb,
                "canon_median": canon,
                "rmse_ratio_vs_1880": rmse_ratio,
                "soft_pass": pareto_soft_pass,
                "verdict": d1884.get("verdict", "UNKNOWN"),
            },
            "qw1885": {
                "robust_score": robust_score,
                "success_rate": succ,
                "rmse_penalty_median": rmse_pen,
                "nonboundary_gain_median": nb_gain,
                "pass": robust_pass,
                "verdict": d1885.get("verdict", "UNKNOWN"),
            },
        },
        "global_score": float(global_score),
        "hard_gate": hard_gate,
        "readiness": readiness,
        "verdict": "NODE_STATE_PHASE7_GATE_COMPLETE",
        "required_next_step": (
            "QW-1887_SIGNED_COUPLING_MICRODYNAMICS_REBUILD" if not robust_pass else "QW-1887_EMPIRICAL_READINESS_BRIDGE"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1886: NODE-STATE PHASE-7 GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        f"- readiness: **{readiness}**",
        f"- hard_gate: **{'PASS' if hard_gate else 'FAIL'}**",
        f"- global_score: {global_score:.3f}",
        "",
        "## Components",
        f"- QW-1882 strict OOS: score={float(s1882):.3f}, pass={strict_pass}",
        f"- QW-1884 pareto: score={pareto_score:.3f}, soft_pass={pareto_soft_pass}",
        f"- QW-1885 robustness: score={robust_score:.3f}, pass={robust_pass}",
        "",
        "## Key Metrics",
        f"- 1884 test: canon={canon:.3f}, nonboundary={nb:.3f}, rmse_ratio_vs_1880={rmse_ratio:.3f}",
        f"- 1885 multisplit: success_rate={succ:.3f}, rmse_penalty_median={rmse_pen:.3f}, nonboundary_gain_median={nb_gain:.3f}",
        "",
        "## Next Required Step",
        f"- {out['required_next_step']}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1886] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1886] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
