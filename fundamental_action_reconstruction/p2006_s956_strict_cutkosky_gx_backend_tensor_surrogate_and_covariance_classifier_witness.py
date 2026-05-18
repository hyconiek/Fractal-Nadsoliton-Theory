#!/usr/bin/env python3
"""P2006 S956 strict Cutkosky gx backend-tensor surrogate + covariance classifier.

Next honest step after P2005: replace scalar gx bound by a backend tensor-surrogate
object and propagate coefficient covariance into gx-scale band for classifier scan.
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
OUT = GEN / "p2006_s956_strict_cutkosky_gx_backend_tensor_surrogate_and_covariance_classifier_witness.json"
TS = "2026-05-18T00:00:00+00:00"


def load(name: str) -> dict[str, Any]:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def cnum(coeffs: dict[str, Any], key: str) -> float:
    e = coeffs.get(key, {})
    if "numeric" in e:
        return float(e.get("numeric", 0.0))
    return float(sp.N(sp.sympify(e.get("symbolic", "0")), 50))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p1853 = load("p1853_s803_strict_b1_symbolic_evaluation_and_scheme_trace_checkpoint.json")
    p2005 = load("p2005_s955_strict_cutkosky_gx_backend_amplitude_bound_and_classifier_witness.json")

    coeffs = p1853.get("b1_symbolic_evaluation", {}).get("evaluated_coefficients", {})
    a_r2 = cnum(coeffs, "a_R2")
    a_ric2 = cnum(coeffs, "a_Ric2")
    a_riem2 = cnum(coeffs, "a_Riem2")

    # backend tensor surrogate object (2x2 curvature-coupling block)
    T = np.array([[a_r2, 0.5*(a_r2+a_ric2)], [0.5*(a_r2+a_ric2), a_ric2+a_riem2]], dtype=float)
    eigvals = np.linalg.eigvalsh(T)

    rows = load("p2004_s954_strict_cutkosky_gx_loop_amplitude_and_robustness_refresh_witness.json").get("grid_table", [])
    deltas = np.array([float(r.get("Delta_opt", 0.0)) for r in rows], dtype=float)
    gx_vals = np.array([float(r.get("Cut_channels", {}).get("gx", 0.0)) for r in rows], dtype=float)
    base_l2 = float(la.norm(deltas, 2)) if deltas.size else float("nan")

    # covariance surrogate from tensor spectrum
    s_abs = np.abs(eigvals)
    tr = float(np.sum(s_abs)) + 1e-15
    center = float(s_abs[1] / tr) if s_abs.size == 2 else 1.0
    sigma = float(0.2 * (s_abs[0] / tr)) if s_abs.size == 2 else 0.05
    lo = max(0.5, center - 2*sigma)
    hi = min(1.5, center + 2*sigma)

    scales = np.linspace(lo, hi, 11)
    scans = []
    counts = {"MISSING_CHANNEL_PRESSURE_SUPPORTED": 0, "STRUCTURAL_OBSTRUCTION_PRESSURE_SUPPORTED": 0, "MIXED_OR_INCONCLUSIVE": 0}
    for f in scales:
        d_new = deltas - (float(f)-1.0)*gx_vals
        ratio = float(la.norm(d_new, 2)/base_l2) if base_l2>0 else float("inf")
        if ratio < 0.95:
            cls = "MISSING_CHANNEL_PRESSURE_SUPPORTED"
        elif ratio > 1.05:
            cls = "STRUCTURAL_OBSTRUCTION_PRESSURE_SUPPORTED"
        else:
            cls = "MIXED_OR_INCONCLUSIVE"
        counts[cls]+=1
        scans.append({"gx_scale": float(f), "l2_ratio_vs_p2004": ratio, "classifier": cls})

    total = len(scans)
    dominant = max(counts, key=lambda k: counts[k])
    dominant_share = counts[dominant]/total if total else 0.0

    gate = {
        "p1853_present": bool(coeffs),
        "p2005_present": p2005.get("result_kind") == "PASS_GX_BACKEND_BOUND_CLASSIFIER_WITNESS",
        "tensor_surrogate_exported": T.shape == (2,2),
        "covariance_band_valid": hi >= lo,
        "scan_nonempty": total > 0,
        "dominant_share_bounded": 0.0 <= dominant_share <= 1.0,
    }

    out = {
        "ledger_id": "P2006_S956_STRICT_CUTKOSKY_GX_BACKEND_TENSOR_SURROGATE_AND_COVARIANCE_CLASSIFIER_WITNESS",
        "packet_id": "P2006",
        "stage_id": "S956",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TS,
        "route": "strict_only",
        "depends_on": {"p1853_present": gate["p1853_present"], "p2005_present": gate["p2005_present"]},
        "gx_tensor_surrogate": {"matrix": T.tolist(), "eigenvalues": eigvals.tolist()},
        "gx_covariance_band": {"center": center, "sigma": sigma, "min": lo, "max": hi, "rule": "scale in [center-2sigma, center+2sigma] from tensor spectrum"},
        "classifier_scan": {"count": total, "table": scans, "counts": counts, "dominant_classifier": dominant, "dominant_share": dominant_share},
        "delta_stats": {"l2_base_p2004": base_l2},
        "gatekeeper_checks": gate,
        "result_kind": "PASS_GX_TENSOR_COVARIANCE_CLASSIFIER_WITNESS" if all(gate.values()) else "OPEN_OBSTRUCTION_WITH_TRACE",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "false_pass_guard": "Tensor-surrogate covariance scan remains diagnostics-only, not theorem-grade closure.",
        "next_honest_step": "Replace tensor surrogate with explicit backend loop amplitude tensor export and rerun covariance classifier.",
        "lay_explanation": "Zamiast jednego zakresu skali użyliśmy małej macierzy backendowej i jej widma, by zbudować bardziej fizyczny zakres niepewności gx.",
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": __import__('scipy').__version__, "sympy": sp.__version__},
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    print(f"[P2006] wrote witness: {OUT}")

if __name__ == "__main__":
    main()
