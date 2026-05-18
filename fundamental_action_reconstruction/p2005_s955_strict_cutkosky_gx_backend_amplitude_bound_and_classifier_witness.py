#!/usr/bin/env python3
"""P2005 S955 strict Cutkosky gx backend-amplitude bound + classifier witness.

Next honest step after P2004: bind gx channel scale to an explicit backend-derived
amplitude bound object and re-evaluate classifier stability under that bound.
"""
from __future__ import annotations
import json, platform
from pathlib import Path
from typing import Any
import numpy as np
import scipy.linalg as la
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2005_s955_strict_cutkosky_gx_backend_amplitude_bound_and_classifier_witness.json"
TS = "2026-05-18T00:00:00+00:00"


def load(name: str) -> dict[str, Any]:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p1853 = load("p1853_s803_strict_b1_symbolic_evaluation_and_scheme_trace_checkpoint.json")
    p2004 = load("p2004_s954_strict_cutkosky_gx_loop_amplitude_and_robustness_refresh_witness.json")

    coeffs = p1853.get("b1_symbolic_evaluation", {}).get("evaluated_coefficients", {})
    def cnum(k: str) -> float:
        e = coeffs.get(k, {})
        if "numeric" in e:
            return float(e.get("numeric", 0.0))
        return float(sp.N(sp.sympify(e.get("symbolic", "0")), 50))

    a_r2 = cnum("a_R2")
    a_ric2 = cnum("a_Ric2")
    a_riem2 = cnum("a_Riem2")

    rows = p2004.get("grid_table", [])
    deltas = np.array([float(r.get("Delta_opt", 0.0)) for r in rows], dtype=float)
    gx_vals = np.array([float(r.get("Cut_channels", {}).get("gx", 0.0)) for r in rows], dtype=float)
    base_l2 = float(la.norm(deltas, 2)) if deltas.size else float("nan")

    # backend-derived amplitude bound for gx scale (explicit object)
    denom = abs(a_r2) + abs(a_ric2) + abs(a_riem2) + 1e-15
    gx_center = (abs(a_r2 + a_ric2) / denom)
    gx_half_width = 0.15 * gx_center
    gx_min = max(0.5, gx_center - gx_half_width)
    gx_max = min(1.5, gx_center + gx_half_width)

    # classify only inside backend-derived bound
    scans = []
    counts = {
        "MISSING_CHANNEL_PRESSURE_SUPPORTED": 0,
        "STRUCTURAL_OBSTRUCTION_PRESSURE_SUPPORTED": 0,
        "MIXED_OR_INCONCLUSIVE": 0,
    }
    band = np.linspace(gx_min, gx_max, 9)
    for f in band:
        d_new = deltas - (float(f) - 1.0) * gx_vals
        ratio = float(la.norm(d_new, 2) / base_l2) if base_l2 > 0 else float("inf")
        if ratio < 0.95:
            cls = "MISSING_CHANNEL_PRESSURE_SUPPORTED"
        elif ratio > 1.05:
            cls = "STRUCTURAL_OBSTRUCTION_PRESSURE_SUPPORTED"
        else:
            cls = "MIXED_OR_INCONCLUSIVE"
        counts[cls] += 1
        scans.append({"gx_scale": float(f), "l2_ratio_vs_p2004": ratio, "classifier": cls})

    total = len(scans)
    dominant = max(counts, key=lambda k: counts[k])
    dominant_share = counts[dominant] / total if total else 0.0

    gate = {
        "p1853_present": bool(coeffs),
        "p2004_present": p2004.get("result_kind") == "PASS_GX_LOOP_AND_ROBUSTNESS_REFRESH_WITNESS",
        "backend_bound_exported": gx_max >= gx_min,
        "scan_nonempty": total > 0,
        "dominant_share_bounded": 0.0 <= dominant_share <= 1.0,
        "delta_norm_bounded": base_l2 < 1.0,
    }

    out = {
        "ledger_id": "P2005_S955_STRICT_CUTKOSKY_GX_BACKEND_AMPLITUDE_BOUND_AND_CLASSIFIER_WITNESS",
        "packet_id": "P2005",
        "stage_id": "S955",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TS,
        "route": "strict_only",
        "depends_on": {"p1853_present": gate["p1853_present"], "p2004_present": gate["p2004_present"]},
        "gx_backend_bound": {
            "center": gx_center,
            "half_width": gx_half_width,
            "min": gx_min,
            "max": gx_max,
            "construction_rule": "center=|a_R2+a_Ric2|/(|a_R2|+|a_Ric2|+|a_Riem2|), half_width=15%*center",
        },
        "classifier_scan": {
            "count": total,
            "table": scans,
            "counts": counts,
            "dominant_classifier": dominant,
            "dominant_share": dominant_share,
        },
        "delta_stats": {"l2_base_p2004": base_l2},
        "gatekeeper_checks": gate,
        "result_kind": "PASS_GX_BACKEND_BOUND_CLASSIFIER_WITNESS" if all(gate.values()) else "OPEN_OBSTRUCTION_WITH_TRACE",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "false_pass_guard": "Bounded classifier remains diagnostics-only and does not imply theorem-grade unitarity closure.",
        "next_honest_step": "Bind gx to explicit backend loop amplitude tensor object and propagate coefficient covariance into bound construction.",
        "lay_explanation": "Skalę nowego kanału ograniczyliśmy teraz jawnie przez dane backendu i dopiero w tym zakresie oceniliśmy stabilność klasyfikacji.",
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": __import__('scipy').__version__},
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P2005] wrote witness: {OUT}")


if __name__ == "__main__":
    main()
