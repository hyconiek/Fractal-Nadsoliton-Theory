#!/usr/bin/env python3
"""P2016 S966 strict Cutkosky channelwise uncertainty transport witness.

Upgrades P2015 by transporting uncertainty at explicit channel level (gg, gh, hh)
instead of a single multiplicative proxy for the entire cut sum.
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
OUT = GEN / "p2016_s966_strict_cutkosky_channelwise_uncertainty_transport_witness.json"
TS = "2026-05-18T00:00:00+00:00"


def load(name: str) -> dict[str, Any]:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p1997 = load("p1997_s947_strict_cutkosky_channelwise_statesum_solver_witness.json")
    rows = p1997.get("grid_table", [])

    e_gg = np.linspace(-0.05, 0.05, 7)
    e_gh = np.linspace(-0.06, 0.06, 7)
    e_hh = np.linspace(-0.08, 0.08, 7)

    table = []
    max_abs_delta = 0.0
    resid_all_positive = True

    for row in rows:
        im = float(row["ImM"])
        ch = row["Cut_channels"]
        gg = float(ch["gg"])
        gh = float(ch["gh"])
        hh = float(ch["hh"])

        deltas = []
        resid_flags = []

        for a in e_gg:
            for b in e_gh:
                for c in e_hh:
                    cut = gg * (1 + float(a)) + gh * (1 + float(b)) + hh * (1 + float(c))
                    d = im - cut
                    deltas.append(d)
                    resid_flags.extend([(gg * (1 + float(a))) > 0, (gh * (1 + float(b))) > 0, (hh * (1 + float(c))) > 0])

        arr = np.array(deltas, dtype=float)
        lo = float(np.min(arr))
        hi = float(np.max(arr))
        max_abs_delta = max(max_abs_delta, float(np.max(np.abs(arr))))
        all_pos = all(resid_flags)
        resid_all_positive = resid_all_positive and all_pos

        table.append({
            "s": row["s"],
            "delta_interval": [lo, hi],
            "delta_center": float(np.mean(arr)),
            "delta_std": float(np.std(arr)),
            "residue_positive_all_samples": all_pos,
        })

    gate = {
        "p1997_present": p1997.get("result_kind") == "PASS_CHANNELWISE_STATESUM_DELTA_OPT_WITNESS",
        "channelwise_uncertainty_table_exported": len(table) == len(rows) and len(rows) > 0,
        "residue_positive_under_channelwise_scan": resid_all_positive,
        "max_abs_delta_channelwise_bounded": max_abs_delta < 0.25,
    }

    out = {
        "ledger_id": "P2016_S966_STRICT_CUTKOSKY_CHANNELWISE_UNCERTAINTY_TRANSPORT_WITNESS",
        "packet_id": "P2016",
        "stage_id": "S966",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TS,
        "route": "strict_only",
        "channel": "graviton->gauge_gauge",
        "depends_on": {"p1997": p1997.get("ledger_id")},
        "scan_parameters": {
            "eps_gg": [-0.05, 0.05, 7],
            "eps_gh": [-0.06, 0.06, 7],
            "eps_hh": [-0.08, 0.08, 7],
            "samples": int(len(e_gg) * len(e_gh) * len(e_hh)),
        },
        "uncertainty_table": table,
        "max_abs_delta_over_scan": max_abs_delta,
        "gatekeeper_checks": gate,
        "result_kind": "PASS_CHANNELWISE_UNCERTAINTY_TRANSPORT_WITNESS" if all(gate.values()) else "OPEN_OBSTRUCTION_WITH_TRACE",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "false_pass_guard": "Channelwise uncertainty transport is still finite-grid diagnostic evidence; it is not a full all-state Cutkosky closure proof.",
        "next_honest_step": "Import explicit loop-derived channel amplitudes from backend integrals and rerun this transport on an extended energy/RG atlas.",
        "toe_progress": "Moves unitarity diagnostics from global proxy perturbation to channel-resolved uncertainty transport, improving strict-lane evidential granularity toward ToE.",
        "lay_explanation": "To dokładniejszy test: zamiast ruszać cały wynik naraz, ruszamy osobno każdy kanał i patrzymy, jak zmienia się błąd całkowity.",
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "sympy": sp.__version__},
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P2016] wrote witness: {OUT}")


if __name__ == "__main__":
    main()
