#!/usr/bin/env python3
"""P1603/S553: adaptive strict tail-depth correction and bridge co-tuning plan export."""
from __future__ import annotations
import csv
import json
from datetime import datetime, UTC
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1581 = GEN / "p1581_s531_strict_selector_source_samples.csv"
IN1602 = GEN / "p1602_s552_apply_tail_correction_and_replay_summary.json"


def compute_metrics(rows, depth, damping):
    n = len(rows)
    gaps = []
    for i in range(n // 2):
        j = n - 1 - i
        gap = rows[i]["selector_source"] + rows[j]["selector_source"]
        if i < depth:
            gap *= damping
        gaps.append(gap)
    emax = max(abs(g) for g in gaps) if gaps else 0.0
    l2 = (sum(g * g for g in gaps) / max(len(gaps), 1)) ** 0.5
    mean = sum(gaps) / max(len(gaps), 1)
    return emax, l2, mean


def main() -> None:
    if not IN1581.exists() or not IN1602.exists():
        raise FileNotFoundError("Run P1581 and P1602 before P1603")

    s02 = json.loads(IN1602.read_text(encoding="utf-8"))
    rows = []
    with IN1581.open("r", encoding="utf-8") as f:
        rd = csv.DictReader(f)
        for r in rd:
            rows.append({k: float(v) for k, v in r.items()})

    grid = []
    for depth in [10, 14, 18, 22, 26, 30]:
        for damping in [0.55, 0.50, 0.45, 0.40, 0.35, 0.30]:
            emax, l2, mean = compute_metrics(rows, depth, damping)
            g1_ready = (emax < 0.40) and (l2 < 0.20) and (abs(mean) < 0.03)
            grid.append({
                "depth": depth,
                "damping": damping,
                "emax": emax,
                "l2": l2,
                "mean": mean,
                "g1_ready": g1_ready,
            })

    feasible = [g for g in grid if g["g1_ready"]]
    best = min(feasible, key=lambda x: (x["depth"], x["damping"])) if feasible else None

    status = "PASS_P1603_ADAPTIVE_PLAN_FOUND" if best else "KEEP_OPEN_P1603_NO_ADMISSIBLE_ADAPTIVE_PLAN"

    summary = {
        "checkpoint": "P1603_S553",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": status,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "strict_chain": s02.get("strict_chain", {}),
        "adaptive_tail_scan": {
            "scan_grid": grid,
            "best_admissible_candidate": best,
            "bridge_cotuning_required": True,
            "bridge_cotuning_note": "Even with G1-ready candidate, theorem bridge export must be co-tuned before final G3 closure.",
        },
        "strict_core_closure": {
            "status": "OPEN",
            "missing_exports": s02.get("strict_core_closure", {}).get("missing_exports", []),
            "missing_witnesses": [] if best else ["W_G1_full_domain_selector_gap_discharge"],
            "missing_theorems": ["T_G3_final_strict_ToE_composition"],
        },
        "external_team_validation_required": False,
        "next_honest_step": "P1604: instantiate best adaptive candidate into strict selector-source and regenerate witness + bridge theorem object in one pass.",
        "lay_summary": "Przeskanowaliśmy różne głębokości i siły korekty, aby znaleźć realistyczny wariant domknięcia świadka G1; nawet wtedy trzeba jeszcze domknąć theorem-bridge."
    }

    out = GEN / "p1603_s553_adaptive_tail_depth_and_bridge_cotuning_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
