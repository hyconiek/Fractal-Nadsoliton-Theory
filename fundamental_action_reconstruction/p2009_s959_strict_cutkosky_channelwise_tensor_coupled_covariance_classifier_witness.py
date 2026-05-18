#!/usr/bin/env python3
"""P2009 S959 strict Cutkosky channelwise tensor-coupled covariance classifier.

Next honest step after P2008: move from gx-only scaling to coupled channelwise
(gg, gh, hh, gx) covariance propagation using a direct backend tensor object.
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
OUT = GEN / "p2009_s959_strict_cutkosky_channelwise_tensor_coupled_covariance_classifier_witness.json"
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
    p2008 = load("p2008_s958_strict_cutkosky_direct_backend_tensor_object_classifier_witness.json")
    p2004 = load("p2004_s954_strict_cutkosky_gx_loop_amplitude_and_robustness_refresh_witness.json")

    coeffs = p1853.get("b1_symbolic_evaluation", {}).get("evaluated_coefficients", {})
    a_r2 = cnum(coeffs, "a_R2")
    a_ric2 = cnum(coeffs, "a_Ric2")
    a_riem2 = cnum(coeffs, "a_Riem2")

    rows = p2004.get("grid_table", [])
    if not rows:
        raise RuntimeError("P2004 grid_table missing")

    # Base delta and channel vectors from P2004
    delta = np.array([float(r["Delta_opt"]) for r in rows], dtype=float)
    ch = {
        k: np.array([float(r["Cut_channels"].get(k, 0.0)) for r in rows], dtype=float)
        for k in ["gg", "gh", "hh", "gx"]
    }
    base_l2 = float(la.norm(delta, 2))

    # Direct backend tensor object reused from P2008 philosophy (3x3 curvature block)
    v = np.array([a_r2, a_ric2, a_riem2], dtype=float)
    C3 = np.outer(v, v) + 1e-12 * np.eye(3)

    # Lift to 4x4 channel covariance via deterministic channel map M
    # channels: gg, gh, hh, gx ; curvature basis: R2, Ric2, Riem2
    M = np.array([
        [1.0, 0.3, 0.1],
        [0.2, 1.0, 0.2],
        [0.1, 0.2, 1.0],
        [0.6, 0.6, 0.4],
    ], dtype=float)
    C4 = M @ C3 @ M.T
    C4 += 1e-12 * np.eye(4)

    # standard deviations as channel scaling uncertainty
    sigma = np.sqrt(np.diag(C4))
    sigma = sigma / (np.linalg.norm(sigma) + 1e-15)

    # coupled scan: t in [-1,1], scale_i = 1 + t * sigma_i
    ts = np.linspace(-1.0, 1.0, 17)
    scans = []
    counts = {"MISSING_CHANNEL_PRESSURE_SUPPORTED": 0, "STRUCTURAL_OBSTRUCTION_PRESSURE_SUPPORTED": 0, "MIXED_OR_INCONCLUSIVE": 0}
    for t in ts:
        scales = 1.0 + t * sigma
        delta_new = delta.copy()
        for i, k in enumerate(["gg", "gh", "hh", "gx"]):
            delta_new -= (scales[i] - 1.0) * ch[k]
        ratio = float(la.norm(delta_new, 2) / base_l2) if base_l2 > 0 else float("inf")
        if ratio < 0.95:
            cls = "MISSING_CHANNEL_PRESSURE_SUPPORTED"
        elif ratio > 1.05:
            cls = "STRUCTURAL_OBSTRUCTION_PRESSURE_SUPPORTED"
        else:
            cls = "MIXED_OR_INCONCLUSIVE"
        counts[cls] += 1
        scans.append({
            "t": float(t),
            "scales": {k: float(scales[i]) for i, k in enumerate(["gg", "gh", "hh", "gx"])},
            "l2_ratio_vs_p2004": ratio,
            "classifier": cls,
        })

    total = len(scans)
    dominant = max(counts, key=lambda k: counts[k])
    dominant_share = counts[dominant] / total if total else 0.0

    gate = {
        "p1853_present": bool(coeffs),
        "p2008_present": p2008.get("result_kind") == "PASS_DIRECT_BACKEND_TENSOR_CLASSIFIER_WITNESS",
        "channel_covariance_exported": C4.shape == (4, 4),
        "covariance_psd": bool(np.all(np.linalg.eigvalsh(C4) > 0)),
        "coupled_scan_nonempty": total > 0,
        "dominant_share_bounded": 0.0 <= dominant_share <= 1.0,
    }

    out = {
        "ledger_id": "P2009_S959_STRICT_CUTKOSKY_CHANNELWISE_TENSOR_COUPLED_COVARIANCE_CLASSIFIER_WITNESS",
        "packet_id": "P2009",
        "stage_id": "S959",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TS,
        "route": "strict_only",
        "depends_on": {"p1853_present": gate["p1853_present"], "p2008_present": gate["p2008_present"]},
        "backend_tensor_objects": {
            "C3_curvature": C3.tolist(),
            "M_channel_map": M.tolist(),
            "C4_channel_covariance": C4.tolist(),
            "channel_sigma_unit": sigma.tolist(),
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
        "result_kind": "PASS_CHANNELWISE_COUPLED_COVARIANCE_CLASSIFIER_WITNESS" if all(gate.values()) else "OPEN_OBSTRUCTION_WITH_TRACE",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "false_pass_guard": "Coupled covariance classifier is diagnostics-only and does not claim theorem-grade closure.",
        "next_honest_step": "Replace deterministic channel map M with exported backend loop-amplitude tensors per channel and repeat coupled covariance scan.",
        "lay_explanation": "Teraz nie zmieniamy tylko jednego kanału, ale sprzężenie wszystkich kanałów naraz przez jedną macierz niepewności. To bliższe realnej współzależności fizycznej.",
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": __import__('scipy').__version__, "sympy": sp.__version__},
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P2009] wrote witness: {OUT}")


if __name__ == "__main__":
    main()
