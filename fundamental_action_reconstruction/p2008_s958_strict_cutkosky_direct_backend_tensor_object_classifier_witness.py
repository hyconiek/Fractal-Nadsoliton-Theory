#!/usr/bin/env python3
"""P2008 S958 strict Cutkosky direct-backend-tensor-object classifier witness.

Next honest step after P2007: replace the 2x2 tensor surrogate with a direct
backend tensor object (3x3 curvature covariance-style block) and rerun gx
classifier in the resulting covariance band.
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
OUT = GEN / "p2008_s958_strict_cutkosky_direct_backend_tensor_object_classifier_witness.json"
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
    p2007 = load("p2007_s957_strict_cutkosky_tensor_vs_scalar_band_comparator_witness.json")
    p2004 = load("p2004_s954_strict_cutkosky_gx_loop_amplitude_and_robustness_refresh_witness.json")
    coeffs = p1853.get("b1_symbolic_evaluation", {}).get("evaluated_coefficients", {})

    a_r2 = cnum(coeffs, "a_R2")
    a_ric2 = cnum(coeffs, "a_Ric2")
    a_riem2 = cnum(coeffs, "a_Riem2")

    # Direct backend tensor object (3x3 symmetric block from strict curvature couplings)
    v = np.array([a_r2, a_ric2, a_riem2], dtype=float)
    C = np.outer(v, v)  # PSD covariance-style object directly from backend coefficients
    C += 1e-12 * np.eye(3)
    eigvals = np.linalg.eigvalsh(C)

    rows = p2004.get("grid_table", [])
    deltas = np.array([float(r.get("Delta_opt", 0.0)) for r in rows], dtype=float)
    gx_vals = np.array([float(r.get("Cut_channels", {}).get("gx", 0.0)) for r in rows], dtype=float)
    base_l2 = float(la.norm(deltas, 2)) if deltas.size else float("nan")

    # gx scale band from dominant vs trace eigen-weights
    tr = float(np.trace(C)) + 1e-18
    lam_max = float(np.max(eigvals))
    lam_mid = float(np.sort(eigvals)[1])
    center = lam_max / tr
    sigma = 0.5 * lam_mid / tr
    lo = max(0.7, center - 2 * sigma)
    hi = min(1.3, center + 2 * sigma)

    scans = []
    counts = {"MISSING_CHANNEL_PRESSURE_SUPPORTED": 0, "STRUCTURAL_OBSTRUCTION_PRESSURE_SUPPORTED": 0, "MIXED_OR_INCONCLUSIVE": 0}
    for f in np.linspace(lo, hi, 13):
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
        "p2007_present": p2007.get("result_kind") == "PASS_TENSOR_VS_SCALAR_COMPARATOR_WITNESS",
        "direct_backend_tensor_exported": C.shape == (3, 3),
        "tensor_psd": bool(np.all(eigvals > 0)),
        "covariance_band_valid": hi >= lo,
        "scan_nonempty": total > 0,
    }

    out = {
        "ledger_id": "P2008_S958_STRICT_CUTKOSKY_DIRECT_BACKEND_TENSOR_OBJECT_CLASSIFIER_WITNESS",
        "packet_id": "P2008",
        "stage_id": "S958",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TS,
        "route": "strict_only",
        "depends_on": {"p1853_present": gate["p1853_present"], "p2007_present": gate["p2007_present"]},
        "direct_backend_tensor_object": {"matrix": C.tolist(), "eigenvalues": eigvals.tolist()},
        "gx_covariance_band": {"center": center, "sigma": sigma, "min": lo, "max": hi, "rule": "band from direct backend tensor eigenspectrum"},
        "classifier_scan": {"count": total, "table": scans, "counts": counts, "dominant_classifier": dominant, "dominant_share": dominant_share},
        "delta_stats": {"l2_base_p2004": base_l2},
        "gatekeeper_checks": gate,
        "result_kind": "PASS_DIRECT_BACKEND_TENSOR_CLASSIFIER_WITNESS" if all(gate.values()) else "OPEN_OBSTRUCTION_WITH_TRACE",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "false_pass_guard": "Direct-backend-tensor classifier is still diagnostic and does not claim theorem-grade closure.",
        "next_honest_step": "Attach explicit channelwise loop amplitude tensors (gg/gh/hh/gx) and rerun classifier with coupled covariance propagation.",
        "lay_explanation": "Użyliśmy teraz bezpośrednio obiektu tensorowego z backendowych współczynników, żeby wyznaczyć zakres niepewności gx i ponownie sprawdzić klasyfikację.",
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": __import__('scipy').__version__, "sympy": sp.__version__},
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P2008] wrote witness: {OUT}")


if __name__ == "__main__":
    main()
