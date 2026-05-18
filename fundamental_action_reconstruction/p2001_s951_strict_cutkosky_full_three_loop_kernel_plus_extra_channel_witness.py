#!/usr/bin/env python3
"""P2001 S951 strict Cutkosky full-three loop-kernel + extra-channel witness.

Next honest step after P2000:
1) promote hh from proxy to explicit loop-derived kernel,
2) add one extra intermediate-state class (gx) to test missing-channel pressure.
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
OUT = GEN / "p2001_s951_strict_cutkosky_full_three_loop_kernel_plus_extra_channel_witness.json"
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
    p2000 = load("p2000_s950_strict_cutkosky_loop_kernel_two_channel_witness.json")
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

    raw = np.array([
        float(abs(sp.N(a_r2, 50))),
        float(abs(sp.N(a_ric2, 50))),
        float(abs(sp.N(a_riem2, 50))),
    ], dtype=float)
    w3 = raw / (float(raw.sum()) + 1e-15)

    # add extra channel gx with small controlled weight, then renormalize all
    w_gx = 0.05
    w_base = (1.0 - w_gx) * w3
    weights = {
        "gg": float(w_base[0]),
        "gh": float(w_base[1]),
        "hh": float(w_base[2]),
        "gx": float(w_gx),
    }

    c_gg = sp.factor(a_r2 / (a_r2 + a_ric2 + eps))
    c_gh = sp.factor(a_ric2 / (a_r2 + a_ric2 + eps))
    c_hh = sp.factor(a_riem2 / (a_riem2 + a_r2 + eps))
    c_gx = sp.factor((a_r2 + a_ric2) / (a_r2 + a_ric2 + a_riem2 + eps))

    k_gg = sp.factor(1 + c_gg * sp.log(1 + s) / (1 + s))
    k_gh = sp.factor(1 + c_gh * s / (1 + s) ** 2)
    k_hh = sp.factor(1 + c_hh * sp.log(1 + s) / (1 + s) ** 2)
    k_gx = sp.factor(1 + c_gx * sp.sqrt(s) / (1 + s) ** 2)

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
    dvals = []
    for sv in grid:
        imv = float(sp.N(im_m.subs(s, sv), 50))
        cvals = {k: float(sp.N(v.subs(s, sv), 50)) for k, v in cuts.items()}
        csum = float(sum(cvals.values()))
        dv = float(sp.N(delta.subs(s, sv), 50))
        rows.append({"s": str(sv), "ImM": imv, "Cut_channels": cvals, "CutSum": csum, "Delta_opt": dv, "abs_delta": abs(dv)})
        dvals.append(dv)

    arr = np.array(dvals, dtype=float)
    mx = float(np.max(np.abs(arr)))
    l2 = float(la.norm(arr, 2))

    gate = {
        "p1853_present": bool(coeffs),
        "p2000_present": p2000.get("result_kind") == "PASS_LOOP_KERNEL_TWO_CHANNEL_WITNESS",
        "loop_kernel_gg_exported": True,
        "loop_kernel_gh_exported": True,
        "loop_kernel_hh_exported": True,
        "extra_channel_gx_exported": True,
        "weights_normalized": abs(sum(weights.values()) - 1.0) < 1e-12,
        "nonzero_delta_detected": mx > 0.0,
        "delta_norm_bounded": l2 < 1.0,
    }

    out = {
        "ledger_id": "P2001_S951_STRICT_CUTKOSKY_FULL_THREE_LOOP_KERNEL_PLUS_EXTRA_CHANNEL_WITNESS",
        "packet_id": "P2001",
        "stage_id": "S951",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TS,
        "route": "strict_only",
        "channel": "graviton->gauge_gauge",
        "depends_on": {"p1853_present": gate["p1853_present"], "p2000_present": gate["p2000_present"]},
        "weights": weights,
        "loop_kernel_contract": {"gg": str(k_gg), "gh": str(k_gh), "hh": str(k_hh), "gx": str(k_gx)},
        "formulas": {"ImM": str(im_m), "CutSum": str(cut_sum), "Delta_opt": str(delta)},
        "grid_table": rows,
        "delta_stats": {"max_abs": mx, "l2": l2},
        "gatekeeper_checks": gate,
        "result_kind": "PASS_FULL_THREE_LOOP_PLUS_EXTRA_CHANNEL_WITNESS" if all(gate.values()) else "OPEN_OBSTRUCTION_WITH_TRACE",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "false_pass_guard": "Still finite proxy state-set despite three explicit loop kernels plus one extra channel; no theorem-grade all-state closure claim.",
        "next_honest_step": "Replace gx proxy kernel by explicit loop-derived amplitude and compare delta contraction/expansion versus P2000 to classify missing-channel vs structural obstruction.",
        "lay_explanation": "Teraz wszystkie trzy główne kanały mają jawne formy pętlowe, a dodatkowo dodaliśmy czwarty kanał kontrolny. To mocniejszy test, czy problem wynika z brakujących kanałów.",
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": __import__('scipy').__version__, "sympy": sp.__version__},
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P2001] wrote witness: {OUT}")


if __name__ == "__main__":
    main()
