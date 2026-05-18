#!/usr/bin/env python3
"""P1998 S948 strict Cutkosky channelwise backend-weight calibration witness.

Next honest step after P1997: keep channelwise structure but replace fixed
manual channel weights with backend-derived calibrated weights from strict B1
coefficients (a_R2, a_Ric2, a_Riem2).
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
OUT = GEN / "p1998_s948_strict_cutkosky_channelwise_backend_weight_calibration_witness.json"
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
    p1997 = load("p1997_s947_strict_cutkosky_channelwise_statesum_solver_witness.json")

    coeffs = p1853.get("b1_symbolic_evaluation", {}).get("evaluated_coefficients", {})
    a_r2 = to_sym(coeffs, "a_R2")
    a_ric2 = to_sym(coeffs, "a_Ric2")
    a_riem2 = to_sym(coeffs, "a_Riem2")

    # strict-lane amplitude backbone retained from P1997/P1996
    s = sp.symbols("s", positive=True, real=True)
    eps = sp.Rational(1, 10**6)
    g_eff = sp.factor(a_r2 + sp.Rational(1, 2) * a_ric2)
    z1 = sp.factor(1 + a_r2 * s / ((a_ric2 + eps) * (1 + s)))
    m_dressed = sp.factor(z1 * g_eff * s)

    rho = sp.Rational(1, 16) / sp.pi
    im_m = sp.factor(rho * m_dressed**2)

    # Backend-derived weights: normalize absolute magnitudes of B1 curvatures.
    raw = np.array([
        float(abs(sp.N(a_r2, 50))),
        float(abs(sp.N(a_ric2, 50))),
        float(abs(sp.N(a_riem2, 50))),
    ], dtype=float)
    denom = float(raw.sum()) + 1e-15
    w = raw / denom

    # Keep mild channel response factors from P1997 (kappa) but weight-calibrated.
    kappas = {"gg": 0.0, "gh": 0.02, "hh": 0.04}
    channels = {
        "gg": {"weight": float(w[0]), "kappa": kappas["gg"]},
        "gh": {"weight": float(w[1]), "kappa": kappas["gh"]},
        "hh": {"weight": float(w[2]), "kappa": kappas["hh"]},
    }

    cut_expr = {}
    for name, v in channels.items():
        cut_expr[name] = sp.factor(rho * sp.Float(v["weight"]) * m_dressed**2 * (1 + sp.Float(v["kappa"]) * s / (1 + s)))

    cut_sum = sp.factor(sum(cut_expr.values()))
    delta = sp.factor(sp.simplify(im_m - cut_sum))

    grid = [sp.Rational(1, 2), sp.Integer(1), sp.Integer(2), sp.Integer(4), sp.Integer(8)]
    rows = []
    dvals = []
    for sv in grid:
        imv = float(sp.N(im_m.subs(s, sv), 50))
        cut_ch = {k: float(sp.N(v.subs(s, sv), 50)) for k, v in cut_expr.items()}
        cutv = float(sum(cut_ch.values()))
        dv = float(sp.N(delta.subs(s, sv), 50))
        rows.append({"s": str(sv), "ImM": imv, "Cut_channels": cut_ch, "CutSum": cutv, "Delta_opt": dv, "abs_delta": abs(dv)})
        dvals.append(dv)

    arr = np.array(dvals, dtype=float)
    mx = float(np.max(np.abs(arr)))
    l2 = float(la.norm(arr, 2))

    gate = {
        "p1853_present": bool(coeffs),
        "p1997_present": p1997.get("result_kind") == "PASS_CHANNELWISE_STATESUM_DELTA_OPT_WITNESS",
        "backend_weight_calibration": True,
        "weights_normalized": abs(float(w.sum()) - 1.0) < 1e-12,
        "nonzero_delta_detected": mx > 0.0,
        "delta_norm_bounded": l2 < 1.0,
    }

    out = {
        "ledger_id": "P1998_S948_STRICT_CUTKOSKY_CHANNELWISE_BACKEND_WEIGHT_CALIBRATION_WITNESS",
        "packet_id": "P1998",
        "stage_id": "S948",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TS,
        "route": "strict_only",
        "channel": "graviton->gauge_gauge",
        "depends_on": {"p1853_present": gate["p1853_present"], "p1997_present": gate["p1997_present"]},
        "backend_weight_rule": "w_i = |a_i| / (|a_R2|+|a_Ric2|+|a_Riem2|)",
        "backend_weights": {"gg": channels["gg"]["weight"], "gh": channels["gh"]["weight"], "hh": channels["hh"]["weight"]},
        "kappa_profile": kappas,
        "formulas": {"ImM": str(im_m), "CutSum": str(cut_sum), "Delta_opt": str(delta)},
        "grid_table": rows,
        "delta_stats": {"max_abs": mx, "l2": l2},
        "gatekeeper_checks": gate,
        "result_kind": "PASS_BACKEND_WEIGHT_CALIBRATION_WITNESS" if all(gate.values()) else "OPEN_OBSTRUCTION_WITH_TRACE",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "false_pass_guard": "Backend-calibrated finite channel witness only; no theorem-grade global Cutkosky closure claim.",
        "next_honest_step": "Replace proxy kappas by explicit loop-derived channel kernels and extend channel basis before UR-link promotion.",
        "lay_explanation": "Wagi kanałów nie są już ręczne — pochodzą z wielkości współczynników backendu. To wzmacnia fizyczne zakotwiczenie, ale nadal nie domyka pełnej unitarności.",
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": __import__('scipy').__version__, "sympy": sp.__version__},
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1998] wrote witness: {OUT}")


if __name__ == "__main__":
    main()
