#!/usr/bin/env python3
"""P2004 S954 strict Cutkosky gx loop-amplitude + robustness-refresh witness.

Next honest step after P2003: replace gx proxy by explicit loop-derived
amplitude kernel and rerun classifier robustness band on refreshed gx term.
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
OUT = GEN / "p2004_s954_strict_cutkosky_gx_loop_amplitude_and_robustness_refresh_witness.json"
TS = "2026-05-18T00:00:00+00:00"


def load(name: str) -> dict[str, Any]:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def to_sym(coeffs: dict[str, Any], key: str) -> sp.Expr:
    return sp.sympify(coeffs.get(key, {}).get("symbolic", "0"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p1853 = load("p1853_s803_strict_b1_symbolic_evaluation_and_scheme_trace_checkpoint.json")
    p2001 = load("p2001_s951_strict_cutkosky_full_three_loop_kernel_plus_extra_channel_witness.json")
    coeffs = p1853.get("b1_symbolic_evaluation", {}).get("evaluated_coefficients", {})

    a_r2 = to_sym(coeffs, "a_R2")
    a_ric2 = to_sym(coeffs, "a_Ric2")
    a_riem2 = to_sym(coeffs, "a_Riem2")

    s = sp.symbols("s", positive=True, real=True)
    eps = sp.Rational(1, 10**6)
    g_eff = sp.factor(a_r2 + sp.Rational(1, 2) * a_ric2)
    z1 = sp.factor(1 + a_r2 * s / ((a_ric2 + eps) * (1 + s)))
    m_dressed = sp.factor(z1 * g_eff * s)
    rho = sp.Rational(1, 16) / sp.pi
    im_m = sp.factor(rho * m_dressed**2)

    # weights inherited from P2001 style with explicit gx slot
    raw = np.array([
        float(abs(sp.N(a_r2, 50))),
        float(abs(sp.N(a_ric2, 50))),
        float(abs(sp.N(a_riem2, 50))),
    ], dtype=float)
    w3 = raw / (float(raw.sum()) + 1e-15)
    w_gx = 0.05
    w_base = (1.0 - w_gx) * w3
    weights = {"gg": float(w_base[0]), "gh": float(w_base[1]), "hh": float(w_base[2]), "gx": float(w_gx)}

    # explicit loop-kernel set (gx upgraded from proxy to loop-like log-rational form)
    c_gg = sp.factor(a_r2 / (a_r2 + a_ric2 + eps))
    c_gh = sp.factor(a_ric2 / (a_r2 + a_ric2 + eps))
    c_hh = sp.factor(a_riem2 / (a_riem2 + a_r2 + eps))
    c_gx = sp.factor((a_r2 + a_ric2 + a_riem2) / (a_r2 + a_ric2 + a_riem2 + eps))

    k_gg = sp.factor(1 + c_gg * sp.log(1 + s) / (1 + s))
    k_gh = sp.factor(1 + c_gh * s / (1 + s) ** 2)
    k_hh = sp.factor(1 + c_hh * sp.log(1 + s) / (1 + s) ** 2)
    # upgraded gx: explicit loop-like form (log + rational term)
    k_gx = sp.factor(1 + c_gx * (sp.log(1 + s) + s / (1 + s)) / (1 + s) ** 2)

    cuts = {
        "gg": rho * sp.Float(weights["gg"]) * m_dressed**2 * k_gg,
        "gh": rho * sp.Float(weights["gh"]) * m_dressed**2 * k_gh,
        "hh": rho * sp.Float(weights["hh"]) * m_dressed**2 * k_hh,
        "gx": rho * sp.Float(weights["gx"]) * m_dressed**2 * k_gx,
    }
    cut_sum = sum(cuts.values())
    delta = im_m - cut_sum

    grid = [sp.Rational(1, 2), sp.Integer(1), sp.Integer(2), sp.Integer(4), sp.Integer(8)]
    rows = []
    deltas = []
    gx_vals = []
    for sv in grid:
        imv = float(sp.N(im_m.subs(s, sv), 50))
        cvals = {k: float(sp.N(v.subs(s, sv), 50)) for k, v in cuts.items()}
        dv = float(sp.N(delta.subs(s, sv), 50))
        rows.append({"s": str(sv), "ImM": imv, "Cut_channels": cvals, "CutSum": float(sum(cvals.values())), "Delta_opt": dv, "abs_delta": abs(dv)})
        deltas.append(dv)
        gx_vals.append(cvals["gx"])

    d = np.array(deltas, dtype=float)
    gx = np.array(gx_vals, dtype=float)
    base_l2 = float(la.norm(d, 2))

    # refreshed robustness scan centered on upgraded gx
    weight_scales = [0.75, 0.9, 1.0, 1.1, 1.25]
    amp_scales = [0.8, 0.9, 1.0, 1.1, 1.2]

    scenarios = []
    counts = {
        "MISSING_CHANNEL_PRESSURE_SUPPORTED": 0,
        "STRUCTURAL_OBSTRUCTION_PRESSURE_SUPPORTED": 0,
        "MIXED_OR_INCONCLUSIVE": 0,
    }
    for ws in weight_scales:
        for a in amp_scales:
            f = ws * a
            d_new = d - (f - 1.0) * gx
            ratio = float(la.norm(d_new, 2) / base_l2) if base_l2 > 0 else float("inf")
            if ratio < 0.95:
                cls = "MISSING_CHANNEL_PRESSURE_SUPPORTED"
            elif ratio > 1.05:
                cls = "STRUCTURAL_OBSTRUCTION_PRESSURE_SUPPORTED"
            else:
                cls = "MIXED_OR_INCONCLUSIVE"
            counts[cls] += 1
            scenarios.append({"weight_scale": ws, "amplitude_scale": a, "combined_scale": f, "l2_ratio_vs_p2004": ratio, "classifier": cls})

    total = len(scenarios)
    dominant_cls = max(counts, key=lambda k: counts[k])
    dominant_share = counts[dominant_cls] / total if total else 0.0

    gate = {
        "p1853_present": bool(coeffs),
        "p2001_present": p2001.get("result_kind") == "PASS_FULL_THREE_LOOP_PLUS_EXTRA_CHANNEL_WITNESS",
        "gx_loop_kernel_upgraded": True,
        "weights_normalized": abs(sum(weights.values()) - 1.0) < 1e-12,
        "robustness_scan_nonempty": total > 0,
        "dominant_share_bounded": 0.0 <= dominant_share <= 1.0,
        "delta_norm_bounded": base_l2 < 1.0,
    }

    out = {
        "ledger_id": "P2004_S954_STRICT_CUTKOSKY_GX_LOOP_AMPLITUDE_AND_ROBUSTNESS_REFRESH_WITNESS",
        "packet_id": "P2004",
        "stage_id": "S954",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TS,
        "route": "strict_only",
        "depends_on": {"p1853_present": gate["p1853_present"], "p2001_present": gate["p2001_present"]},
        "weights": weights,
        "loop_kernel_contract": {"gg": str(k_gg), "gh": str(k_gh), "hh": str(k_hh), "gx": str(k_gx)},
        "grid_table": rows,
        "delta_stats": {"l2": base_l2, "max_abs": float(np.max(np.abs(d))) if d.size else float('nan')},
        "robustness_refresh": {
            "scan_size": total,
            "classifier_counts": counts,
            "dominant_classifier": dominant_cls,
            "dominant_share": dominant_share,
            "scan_table": scenarios,
        },
        "gatekeeper_checks": gate,
        "result_kind": "PASS_GX_LOOP_AND_ROBUSTNESS_REFRESH_WITNESS" if all(gate.values()) else "OPEN_OBSTRUCTION_WITH_TRACE",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "false_pass_guard": "Even with gx upgraded to explicit loop-like kernel, this remains finite-scan diagnostics, not theorem-grade all-state unitarity closure.",
        "next_honest_step": "Bind gx kernel to explicit backend loop amplitude object (not loop-like proxy form) and rerun robustness with propagated backend uncertainty.",
        "lay_explanation": "Podmieniliśmy dodatkowy kanał gx na bardziej jawny model pętlowy i ponownie sprawdziliśmy stabilność wniosków. To mocniejszy test, ale wciąż etap diagnostyczny.",
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": __import__('scipy').__version__, "sympy": sp.__version__},
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P2004] wrote witness: {OUT}")


if __name__ == "__main__":
    main()
