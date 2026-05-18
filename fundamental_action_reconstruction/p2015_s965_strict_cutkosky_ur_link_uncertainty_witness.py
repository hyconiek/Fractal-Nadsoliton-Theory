#!/usr/bin/env python3
"""P2015 S965 strict Cutkosky UR-link uncertainty witness.

Constructs a bounded uncertainty table for the strict-lane optical theorem proxy
in graviton -> gauge_gauge channel and checks residue positivity stability under
admissible RG-safe perturbations around the current strict diagnostic model.
"""
from __future__ import annotations

import json
import platform
from pathlib import Path
from typing import Any

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2015_s965_strict_cutkosky_ur_link_uncertainty_witness.json"
TS = "2026-05-18T00:00:00+00:00"


def load(name: str) -> dict[str, Any]:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p1997 = load("p1997_s947_strict_cutkosky_channelwise_statesum_solver_witness.json")

    coeffs = p1997.get("formulas", {})
    if not coeffs:
        payload = {
            "ledger_id": "P2015_S965_STRICT_CUTKOSKY_UR_LINK_UNCERTAINTY_WITNESS",
            "packet_id": "P2015",
            "stage_id": "S965",
            "produced_by": Path(__file__).name,
            "timestamp_utc": TS,
            "result_kind": "OPEN_MISSING_ARTIFACT",
            "status": "OPEN_OBSTRUCTION_WITH_TRACE",
            "false_pass_guard": "Requires upstream artifacts; no closure claim.",
        }
        OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(f"[P2015] wrote witness: {OUT}")
        return

    rows = p1997.get("grid_table", [])

    # RG-safe local perturbation family (diagnostic only):
    # Z1 -> (1+e1) Z1, channel dressing -> (1+e2 s/(1+s)).
    e1_vals = np.linspace(-0.03, 0.03, 7)
    e2_vals = np.linspace(-0.05, 0.05, 11)
    poles = [0.5, 1.0, 2.0]

    table = []
    residue_all_positive = True
    max_abs_delta = 0.0

    for row in rows:
        sv = float(sp.Rational(row["s"]))
        base_im = float(row["ImM"])
        deltas = []
        resid_flags = []

        for e1 in e1_vals:
            for e2 in e2_vals:
                scale = (1.0 + float(e1)) ** 2 * (1.0 + float(e2) * (sv / (1.0 + sv)))
                cut = base_im * scale
                delta = base_im - cut
                deltas.append(delta)

                for p in poles:
                    # proxy residue monitor under multiplicative perturbation
                    rp = base_im * (1.0 + float(e1)) ** 2 * (1.0 + float(e2) * (p / (1.0 + p)))
                    ok = rp > 0
                    resid_flags.append(ok)

        d_arr = np.array(deltas, dtype=float)
        lo = float(np.min(d_arr))
        hi = float(np.max(d_arr))
        max_abs_delta = max(max_abs_delta, float(np.max(np.abs(d_arr))))
        all_pos = all(resid_flags)
        residue_all_positive = residue_all_positive and all_pos

        table.append(
            {
                "s": str(sv),
                "delta_interval": [lo, hi],
                "delta_center": float(np.mean(d_arr)),
                "delta_std": float(np.std(d_arr)),
                "residue_positive_all_samples": all_pos,
            }
        )

    gate = {
        "p1997_present": p1997.get("result_kind") == "PASS_CHANNELWISE_STATESUM_DELTA_OPT_WITNESS",
        "bounded_uncertainty_table_exported": len(table) == len(rows) and len(rows) > 0,
        "residue_positive_under_scan": residue_all_positive,
        "max_abs_delta_small_scan": max_abs_delta < 0.2,
    }

    out = {
        "ledger_id": "P2015_S965_STRICT_CUTKOSKY_UR_LINK_UNCERTAINTY_WITNESS",
        "packet_id": "P2015",
        "stage_id": "S965",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TS,
        "route": "strict_only",
        "channel": "graviton->gauge_gauge",
        "depends_on": {"p1997": p1997.get("ledger_id")},
        "scan_parameters": {
            "epsilon_z1_min": float(e1_vals.min()),
            "epsilon_z1_max": float(e1_vals.max()),
            "epsilon_ch_min": float(e2_vals.min()),
            "epsilon_ch_max": float(e2_vals.max()),
            "samples": int(len(e1_vals) * len(e2_vals)),
        },
        "uncertainty_table": table,
        "max_abs_delta_over_scan": max_abs_delta,
        "gatekeeper_checks": gate,
        "result_kind": "PASS_UR_LINK_UNCERTAINTY_WITNESS" if all(gate.values()) else "OPEN_OBSTRUCTION_WITH_TRACE",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "false_pass_guard": "This is a bounded local uncertainty scan around current strict diagnostic proxies, not a full all-state Cutkosky theorem.",
        "next_honest_step": "Replace proxy multiplicative perturbation with explicit loop-derived channel amplitudes and re-run positivity/uncertainty transport across a wider RG-safe atlas.",
        "toe_progress": "Adds machine-checkable uncertainty and positivity stability evidence for UR-link diagnostics, narrowing one obstruction in strict-lane unitarity work without claiming ToE closure.",
        "lay_explanation": "Sprawdziliśmy, czy przy małych, kontrolowanych zmianach modelu wyniki nie 'psują się' gwałtownie: to zwiększa zaufanie do stabilności, ale nie jest jeszcze pełnym dowodem unitarności.",
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "sympy": sp.__version__},
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P2015] wrote witness: {OUT}")


if __name__ == "__main__":
    main()
