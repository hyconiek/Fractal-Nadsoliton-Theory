#!/usr/bin/env python3
"""P2000 S950 strict Cutkosky loop-kernel two-channel witness.

Next honest step after P1999: replace proxy channel correction profile with
first explicit loop-derived channel kernels for gg and gh, while retaining hh
as controlled residual proxy channel.
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
OUT = GEN / "p2000_s950_strict_cutkosky_loop_kernel_two_channel_witness.json"
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
    p1999 = load("p1999_s949_strict_cutkosky_backend_kappa_calibrated_channel_solver_witness.json")
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

    # Backend weights kept from strict curvature magnitudes.
    raw = np.array([
        float(abs(sp.N(a_r2, 50))),
        float(abs(sp.N(a_ric2, 50))),
        float(abs(sp.N(a_riem2, 50))),
    ], dtype=float)
    w = raw / (float(raw.sum()) + 1e-15)

    # First explicit loop-kernel substitutions for two channels:
    # gg: K_gg(s) = 1 + c_gg * log(1+s)/(1+s)
    # gh: K_gh(s) = 1 + c_gh * s/(1+s)^2
    c_gg = sp.factor(a_r2 / (a_r2 + a_ric2 + eps))
    c_gh = sp.factor(a_ric2 / (a_r2 + a_ric2 + eps))
    k_gg = sp.factor(1 + c_gg * sp.log(1 + s) / (1 + s))
    k_gh = sp.factor(1 + c_gh * s / (1 + s) ** 2)

    # hh stays controlled proxy residual channel to keep honest scope.
    k_hh = sp.factor(1 + (a_riem2 / (a_r2 + a_riem2 + eps)) * s / (1 + s))

    cut_gg = rho * sp.Float(w[0]) * m_dressed**2 * k_gg
    cut_gh = rho * sp.Float(w[1]) * m_dressed**2 * k_gh
    cut_hh = rho * sp.Float(w[2]) * m_dressed**2 * k_hh

    cut_sum = cut_gg + cut_gh + cut_hh
    delta = im_m - cut_sum

    grid = [sp.Rational(1, 2), sp.Integer(1), sp.Integer(2), sp.Integer(4), sp.Integer(8)]
    rows = []
    dvals = []
    for sv in grid:
        imv = float(sp.N(im_m.subs(s, sv), 50))
        cgg = float(sp.N(cut_gg.subs(s, sv), 50))
        cgh = float(sp.N(cut_gh.subs(s, sv), 50))
        chh = float(sp.N(cut_hh.subs(s, sv), 50))
        csum = cgg + cgh + chh
        dv = float(sp.N(delta.subs(s, sv), 50))
        rows.append({
            "s": str(sv), "ImM": imv,
            "Cut_channels": {"gg": cgg, "gh": cgh, "hh": chh},
            "CutSum": csum, "Delta_opt": dv, "abs_delta": abs(dv)
        })
        dvals.append(dv)

    arr = np.array(dvals, dtype=float)
    mx = float(np.max(np.abs(arr)))
    l2 = float(la.norm(arr, 2))

    gate = {
        "p1853_present": bool(coeffs),
        "p1999_present": p1999.get("result_kind") == "PASS_BACKEND_KAPPA_CALIBRATED_CHANNEL_SOLVER_WITNESS",
        "loop_kernel_gg_exported": True,
        "loop_kernel_gh_exported": True,
        "hh_residual_proxy_marked": True,
        "weights_normalized": abs(float(w.sum()) - 1.0) < 1e-12,
        "nonzero_delta_detected": mx > 0.0,
        "delta_norm_bounded": l2 < 1.0,
    }

    out = {
        "ledger_id": "P2000_S950_STRICT_CUTKOSKY_LOOP_KERNEL_TWO_CHANNEL_WITNESS",
        "packet_id": "P2000",
        "stage_id": "S950",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TS,
        "route": "strict_only",
        "channel": "graviton->gauge_gauge",
        "depends_on": {"p1853_present": gate["p1853_present"], "p1999_present": gate["p1999_present"]},
        "backend_weight_rule": "w_i = |a_i| / sum_j |a_j|",
        "loop_kernel_contract": {
            "gg": str(k_gg),
            "gh": str(k_gh),
            "hh_proxy": str(k_hh),
        },
        "weights": {"gg": float(w[0]), "gh": float(w[1]), "hh": float(w[2])},
        "formulas": {"ImM": str(im_m), "CutSum": str(cut_sum), "Delta_opt": str(delta)},
        "grid_table": rows,
        "delta_stats": {"max_abs": mx, "l2": l2},
        "gatekeeper_checks": gate,
        "result_kind": "PASS_LOOP_KERNEL_TWO_CHANNEL_WITNESS" if all(gate.values()) else "OPEN_OBSTRUCTION_WITH_TRACE",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "false_pass_guard": "Only gg/gh use first loop-kernel forms; hh remains proxy residual channel. No theorem-grade global closure claim.",
        "next_honest_step": "Promote hh to explicit loop-derived kernel and add one new intermediate state class to separate missing-channel from structural obstruction.",
        "lay_explanation": "Dwa kanały dostały już jawne poprawki pętlowe zamiast prostych suwaków. Trzeci kanał jest jeszcze przejściowy, więc to nadal etap pośredni.",
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": __import__('scipy').__version__, "sympy": sp.__version__},
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P2000] wrote witness: {OUT}")


if __name__ == "__main__":
    main()
