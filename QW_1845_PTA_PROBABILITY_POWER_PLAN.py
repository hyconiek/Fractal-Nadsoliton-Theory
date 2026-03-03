#!/usr/bin/env python3
"""
QW-1845: PTA probability-threshold power plan.

Designs required replication count for one-sided exact binomial test:
H0: p <= p0 vs H1: p > p0, where p0=0.9 from QW-1839.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Optional

from scipy.stats import binom


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1845_pta_probability_power_plan.json"
OUT_MD = ROOT / "RAPORT_QW1845_PTA_PROBABILITY_POWER_PLAN.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def binom_tail(n: int, k: int, p: float) -> float:
    return float(binom.sf(k - 1, n, p))


def find_critical_k(n: int, p0: float, alpha: float) -> Optional[int]:
    # minimal k such that tail under p0 is <= alpha
    for k in range(0, n + 1):
        if binom_tail(n, k, p0) <= alpha:
            return k
    return None


def min_n_for_power(p0: float, p1: float, alpha: float, power_target: float, n_max: int = 500) -> Optional[Dict[str, float]]:
    for n in range(1, n_max + 1):
        kcrit = find_critical_k(n, p0=p0, alpha=alpha)
        if kcrit is None:
            continue
        power = binom_tail(n, kcrit, p1)
        if power >= power_target:
            return {
                "n": n,
                "k_critical": kcrit,
                "alpha_actual": binom_tail(n, kcrit, p0),
                "power": power,
            }
    return None


def main() -> None:
    d1839 = load("report_qw1839_joint_confirmatory_prereg_protocol.json")
    d1843 = load("report_qw1843_pta_threshold_inference_rigor.json")

    p0 = float(d1839["protocol"]["pta_protocol"]["thresholds"]["prob_quantile_gain_positive_min"])
    n_current = int(d1843.get("n_replications", 0))

    alpha = 0.05
    power_target = 0.80
    p1_grid = [0.93, 0.95, 0.97, 0.99, 1.00]

    rows = []
    for p1 in p1_grid:
        res = min_n_for_power(p0=p0, p1=p1, alpha=alpha, power_target=power_target, n_max=800)
        if res is None:
            rows.append(
                {
                    "p_true_assumed": p1,
                    "feasible": False,
                }
            )
        else:
            n_need = int(res["n"])
            rows.append(
                {
                    "p_true_assumed": p1,
                    "feasible": True,
                    "n_required": n_need,
                    "k_critical": int(res["k_critical"]),
                    "alpha_actual": float(res["alpha_actual"]),
                    "power": float(res["power"]),
                    "additional_needed_vs_current": max(0, n_need - n_current),
                }
            )

    # Conservative recommendation = max required among feasible scenarios up to p_true=0.97
    feas_upto_097 = [r for r in rows if r.get("feasible") and r["p_true_assumed"] <= 0.97]
    if feas_upto_097:
        n_reco = max(int(r["n_required"]) for r in feas_upto_097)
    else:
        n_reco = None

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "test_design": {
            "null": f"p <= {p0}",
            "alternative": f"p > {p0}",
            "alpha": alpha,
            "power_target": power_target,
            "current_n": n_current,
        },
        "grid_results": rows,
        "recommendation": {
            "recommended_total_n_for_p_true_up_to_0p97": n_reco,
            "additional_needed_vs_current": (max(0, n_reco - n_current) if n_reco is not None else None),
        },
        "verdict": "PTA_PROBABILITY_POWER_PLAN_COMPLETE",
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1845: PTA PROBABILITY POWER PLAN",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Test: H0 p<={p0:.3f} vs H1 p>{p0:.3f}",
        f"- alpha={alpha:.3f}, power_target={power_target:.2f}",
        f"- current n={n_current}",
        "",
        "## Grid (p_true assumptions)",
    ]

    for r in rows:
        if not r.get("feasible"):
            lines.append(f"- p_true={r['p_true_assumed']:.2f}: no solution in n<=800")
        else:
            lines.append(
                (
                    f"- p_true={r['p_true_assumed']:.2f}: n={r['n_required']}, "
                    f"k_crit={r['k_critical']}, alpha_act={r['alpha_actual']:.4f}, "
                    f"power={r['power']:.4f}, add_vs_current={r['additional_needed_vs_current']}"
                )
            )

    lines += [
        "",
        "## Recommendation",
        f"- recommended_total_n_for_p_true_up_to_0.97: {n_reco}",
        (
            "- additional_needed_vs_current: "
            f"{out['recommendation']['additional_needed_vs_current']}"
        ),
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1845] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1845] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
