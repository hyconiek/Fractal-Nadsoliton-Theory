#!/usr/bin/env python3
"""
QW-1846: Feasibility map for PTA probability threshold criterion.

Maps minimal sample sizes for one-sided exact binomial testing
across threshold p0 and assumed true p.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Optional

from scipy.stats import binom


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1846_pta_prob_threshold_feasibility_map.json"
OUT_MD = ROOT / "RAPORT_QW1846_PTA_PROB_THRESHOLD_FEASIBILITY_MAP.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def binom_tail(n: int, k: int, p: float) -> float:
    # Survival function gives P[X >= k] as sf(k-1).
    return float(binom.sf(k - 1, n, p))


def find_critical_k(n: int, p0: float, alpha: float) -> Optional[int]:
    for k in range(0, n + 1):
        if binom_tail(n, k, p0) <= alpha:
            return k
    return None


def min_n_for_power(p0: float, p1: float, alpha: float, power_target: float, n_max: int = 400) -> Optional[Dict[str, float]]:
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


def max_p0_rejectable_all_positive(n: int, alpha: float) -> float:
    # for k=n, reject iff p0^n <= alpha
    return float(alpha ** (1.0 / n))


def main() -> None:
    d1843 = load("report_qw1843_pta_threshold_inference_rigor.json")

    n_current = int(d1843.get("n_replications", 0))
    alpha = 0.05
    power_target = 0.80

    p0_grid = [0.75, 0.80, 0.85, 0.90, 0.95]
    p1_grid = [0.90, 0.93, 0.95, 0.97, 0.99, 1.00]

    table = []
    for p0 in p0_grid:
        for p1 in p1_grid:
            if p1 <= p0:
                table.append(
                    {
                        "p0_threshold": p0,
                        "p_true_assumed": p1,
                        "feasible": False,
                        "reason": "p_true_not_above_threshold",
                    }
                )
                continue

            res = min_n_for_power(p0=p0, p1=p1, alpha=alpha, power_target=power_target, n_max=800)
            if res is None:
                table.append(
                    {
                        "p0_threshold": p0,
                        "p_true_assumed": p1,
                        "feasible": False,
                        "reason": "no_solution_in_range",
                    }
                )
            else:
                table.append(
                    {
                        "p0_threshold": p0,
                        "p_true_assumed": p1,
                        "feasible": True,
                        "n_required": int(res["n"]),
                        "k_critical": int(res["k_critical"]),
                        "alpha_actual": float(res["alpha_actual"]),
                        "power": float(res["power"]),
                        "add_vs_current": max(0, int(res["n"]) - n_current),
                    }
                )

    p0_max_at_current_n = max_p0_rejectable_all_positive(n_current, alpha)

    # Highlight row for current criterion p0=0.90, p_true=0.97 as a practical scenario.
    focus = [r for r in table if r.get("p0_threshold") == 0.90 and r.get("p_true_assumed") == 0.97 and r.get("feasible")]
    focus_row = focus[0] if focus else None

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "settings": {
            "alpha": alpha,
            "power_target": power_target,
            "n_current": n_current,
            "p0_grid": p0_grid,
            "p1_grid": p1_grid,
        },
        "current_n_limit": {
            "max_threshold_rejectable_if_all_positive": p0_max_at_current_n,
        },
        "feasibility_table": table,
        "focus_current_criterion_p0_0p90_p1_0p97": focus_row,
        "verdict": "PTA_PROB_THRESHOLD_FEASIBILITY_MAP_COMPLETE",
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1846: PTA PROB THRESHOLD FEASIBILITY MAP",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- alpha={alpha:.3f}, power_target={power_target:.2f}, n_current={n_current}",
        (
            "- Max p0 rejectable at current n (all-positive case): "
            f"{p0_max_at_current_n:.4f}"
        ),
        "",
        "## Focus (current criterion)",
    ]

    if focus_row is None:
        lines.append("- p0=0.90, p_true=0.97: no feasible solution in scanned range")
    else:
        lines.append(
            (
                f"- p0=0.90, p_true=0.97 -> n_required={focus_row['n_required']}, "
                f"k_crit={focus_row['k_critical']}, add_vs_current={focus_row['add_vs_current']}"
            )
        )

    lines += [
        "",
        "## Selected rows",
    ]

    for r in table:
        if r.get("p0_threshold") in (0.80, 0.85, 0.90) and r.get("p_true_assumed") in (0.95, 0.97, 0.99):
            if r.get("feasible"):
                lines.append(
                    (
                        f"- p0={r['p0_threshold']:.2f}, p_true={r['p_true_assumed']:.2f}: "
                        f"n={r['n_required']}, add={r['add_vs_current']}"
                    )
                )
            else:
                lines.append(
                    (
                        f"- p0={r['p0_threshold']:.2f}, p_true={r['p_true_assumed']:.2f}: "
                        f"not feasible ({r.get('reason')})"
                    )
                )

    lines += [
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1846] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1846] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
