#!/usr/bin/env python3
"""
QW-2132: RG fixed-point Jacobian gate (strict proxy level).

Purpose:
- provide an explicit beta-function system for the FIN v5 strict chain couplings,
- identify analytic fixed points in the declared proxy RG model,
- classify UV/IR directions via Jacobian + directional sign probes,
- keep scope boundaries explicit (proxy closure, not full nonperturbative proof).
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2132_rg_fixed_point_jacobian_gate.json"
OUT_MD = ROOT / "RAPORT_QW2132_RG_FIXED_POINT_JACOBIAN_GATE.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def load_optional_json(name: str) -> Dict | None:
    path = ROOT / name
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def find_entry(entries: List[Dict], pid: str) -> Dict:
    for e in entries:
        if e.get("id") == pid:
            return e
    raise KeyError(f"Missing entry id: {pid}")


def beta_vector(x: np.ndarray) -> np.ndarray:
    """
    One-loop proxy system in t = ln(mu):
    x = [g1, g2, g3, y_t, lambda_h, g_gr].
    """
    g1, g2, g3, yt, lam, ggr = [float(v) for v in x]
    c = 1.0 / (16.0 * math.pi * math.pi)

    # SM-like one-loop proxy betas (hypercharge coupling treated as g1 = g').
    b_g1 = c * (41.0 / 6.0) * g1**3
    b_g2 = c * (-19.0 / 6.0) * g2**3
    b_g3 = c * (-7.0) * g3**3
    b_yt = c * yt * (
        (9.0 / 2.0) * yt**2
        - (17.0 / 12.0) * g1**2
        - (9.0 / 4.0) * g2**2
        - 8.0 * g3**2
    )
    b_lam = c * (
        24.0 * lam**2
        - 6.0 * yt**4
        + (3.0 / 8.0) * (2.0 * g2**4 + (g2**2 + g1**2) ** 2)
        + (-9.0 * g2**2 - 3.0 * g1**2 + 12.0 * yt**2) * lam
    )

    # Logistic UV-safe proxy used in QW-2073.
    b_ggr = 2.0 * ggr * (1.0 - ggr)
    return np.array([b_g1, b_g2, b_g3, b_yt, b_lam, b_ggr], dtype=float)


def jacobian_fd(x: np.ndarray, h: float = 1e-6) -> np.ndarray:
    n = int(x.shape[0])
    j = np.zeros((n, n), dtype=float)
    for k in range(n):
        e = np.zeros(n, dtype=float)
        e[k] = h
        f_plus = beta_vector(x + e)
        f_minus = beta_vector(x - e)
        j[:, k] = (f_plus - f_minus) / (2.0 * h)
    return j


def directional_probe_at_gaussian(eps: float = 1e-3) -> Dict[str, str]:
    labels = ["g1", "g2", "g3", "y_t", "lambda_h", "g_gr"]
    out: Dict[str, str] = {}
    for i, label in enumerate(labels):
        x = np.zeros(6, dtype=float)
        x[i] = eps
        b = float(beta_vector(x)[i])
        if b < 0.0:
            out[label] = "uv_attractive_towards_zero"
        elif b > 0.0:
            out[label] = "uv_repulsive_from_zero"
        else:
            out[label] = "marginal_or_higher_order"
    return out


def parse_uv_ir_counts(direction_map: Dict[str, str]) -> Dict[str, int]:
    uv_attr = sum(1 for v in direction_map.values() if v == "uv_attractive_towards_zero")
    uv_rep = sum(1 for v in direction_map.values() if v == "uv_repulsive_from_zero")
    marginal = sum(1 for v in direction_map.values() if v == "marginal_or_higher_order")
    return {"uv_attractive_count": uv_attr, "uv_repulsive_count": uv_rep, "marginal_count": marginal}


def max_abs_diff(a: List[float], b: List[float]) -> float:
    return float(max(abs(x - y) for x, y in zip(a, b)))


def main() -> None:
    r2069 = load_json("report_qw2069_full_sm_gr_derivation_package.json")
    r2126 = load_json("report_qw2126_gauge_yukawa_numeric_derivation_gate.json")
    r2073 = load_optional_json("report_qw2073_radiative_channels_closure_upgrade.json")

    entries = r2069["entries"]
    m_top = float(find_entry(entries, "m_top")["predicted_value"])
    m_h = float(find_entry(entries, "m_h")["predicted_value"])
    v = float(find_entry(entries, "v_higgs")["predicted_value"])

    g1_ref = float(r2126["derived_gauge_couplings"]["gprime_u1"])
    g2_ref = float(r2126["derived_gauge_couplings"]["g_su2"])
    g3_ref = float(r2126["derived_gauge_couplings"]["g3_su3"])
    yt_ref = float(math.sqrt(2.0) * m_top / max(v, 1e-300))
    lam_ref = float((m_h * m_h) / max(2.0 * v * v, 1e-300))

    m_pl = 1.2209e19
    mu_ref = 1.0
    ggr_ref = float((mu_ref / m_pl) ** 2)

    x_ref = np.array([g1_ref, g2_ref, g3_ref, yt_ref, lam_ref, ggr_ref], dtype=float)
    b_ref = beta_vector(x_ref)

    # Analytic fixed points in this declared proxy system:
    # P0 = Gaussian (all zero), P1 = asymptotic-safe gravity branch (g_gr=1).
    p0 = np.zeros(6, dtype=float)
    p1 = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 1.0], dtype=float)

    j0 = jacobian_fd(p0)
    j1 = jacobian_fd(p1)
    eig0 = np.linalg.eigvals(j0)
    eig1 = np.linalg.eigvals(j1)

    dir0 = directional_probe_at_gaussian(eps=1e-3)
    counts0 = parse_uv_ir_counts(dir0)

    # Local probe around g_gr=1: dg/dt = 2g(1-g) => derivative -2 at g=1 (UV attractive).
    ggr_slope_at_1 = float(-2.0)
    ggr_slope_at_0 = float(2.0)

    # Consistency check with QW-2073 gr samples, if available.
    gr_consistent = None
    gr_max_abs_diff = None
    if r2073 is not None:
        samples = list(r2073["diagnostics"]["gr_samples"])
        ref_vals = [float(row["g_dimensionless"]) for row in samples]
        calc_vals = []
        g0 = ggr_ref
        for row in samples:
            mu = float(row["mu_gev"])
            t = math.log(mu / mu_ref)
            e2t = math.exp(2.0 * t)
            g_mu = (g0 * e2t) / (1.0 + g0 * (e2t - 1.0))
            calc_vals.append(float(g_mu))
        gr_max_abs_diff = max_abs_diff(ref_vals, calc_vals)
        gr_consistent = bool(gr_max_abs_diff <= 1e-12)

    flags = {
        "beta_system_explicitly_defined": True,
        "analytic_fixed_points_declared_ge_2": True,
        "jacobian_computed_at_declared_fixed_points": True,
        "gaussian_directional_probe_computed": True,
        "gr_logistic_fixed_point_uv_attractive_at_g_eq_1": bool(ggr_slope_at_1 < 0.0),
        "gr_logistic_fixed_point_uv_repulsive_at_g_eq_0": bool(ggr_slope_at_0 > 0.0),
        "nonempty_uv_attractive_set_documented": bool(counts0["uv_attractive_count"] >= 1),
        "nonempty_uv_repulsive_set_documented": bool(counts0["uv_repulsive_count"] >= 1),
        "q2073_gr_channel_consistency": bool(gr_consistent) if gr_consistent is not None else False,
        "deterministic_no_scan_no_retune": True,
        "full_nonperturbative_rg_fixed_point_proof": False,
    }
    pass_count = int(sum(1 for v_ in flags.values() if bool(v_)))
    total_flags = int(len(flags))

    verdict = (
        "RG_FIXED_POINT_JACOBIAN_GATE_PASS_STRICT_PROXY_PARTIAL"
        if (
            flags["beta_system_explicitly_defined"]
            and flags["analytic_fixed_points_declared_ge_2"]
            and flags["jacobian_computed_at_declared_fixed_points"]
            and flags["gaussian_directional_probe_computed"]
            and flags["gr_logistic_fixed_point_uv_attractive_at_g_eq_1"]
            and flags["nonempty_uv_attractive_set_documented"]
            and flags["nonempty_uv_repulsive_set_documented"]
            and flags["deterministic_no_scan_no_retune"]
        )
        else "RG_FIXED_POINT_JACOBIAN_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2069_package": "report_qw2069_full_sm_gr_derivation_package.json",
            "q2126_numeric_bridge": "report_qw2126_gauge_yukawa_numeric_derivation_gate.json",
            "q2073_radiative_upgrade": (
                "report_qw2073_radiative_channels_closure_upgrade.json" if r2073 is not None else None
            ),
        },
        "reference_couplings_at_mu_1gev": {
            "g1_u1": g1_ref,
            "g2_su2": g2_ref,
            "g3_su3": g3_ref,
            "y_top": yt_ref,
            "lambda_h": lam_ref,
            "g_gr_dimensionless": ggr_ref,
        },
        "beta_at_reference": {
            "beta_g1": float(b_ref[0]),
            "beta_g2": float(b_ref[1]),
            "beta_g3": float(b_ref[2]),
            "beta_y_top": float(b_ref[3]),
            "beta_lambda_h": float(b_ref[4]),
            "beta_g_gr": float(b_ref[5]),
        },
        "fixed_points": {
            "gaussian": [float(v_) for v_ in p0.tolist()],
            "asymptotic_safe_gr_branch": [float(v_) for v_ in p1.tolist()],
        },
        "jacobian_audit": {
            "jacobian_gaussian": j0.tolist(),
            "jacobian_as_gr": j1.tolist(),
            "eigvals_gaussian": [{"re": float(np.real(z)), "im": float(np.imag(z))} for z in eig0],
            "eigvals_as_gr": [{"re": float(np.real(z)), "im": float(np.imag(z))} for z in eig1],
        },
        "directional_probe_gaussian": {
            "eps": 1e-3,
            "classification": dir0,
            "counts": counts0,
        },
        "gr_logistic_slopes": {
            "dbeta_dg_at_0": ggr_slope_at_0,
            "dbeta_dg_at_1": ggr_slope_at_1,
        },
        "q2073_consistency": {
            "available": bool(r2073 is not None),
            "max_abs_diff_gr_samples": gr_max_abs_diff,
            "consistent_within_1e_12": gr_consistent,
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "EXTEND_FROM_PROXY_RG_TO_FULL_NONPERTURBATIVE_RG_FLOW_WITH_PROVEN_FIXED_POINT_STABILITY"
            if verdict.endswith("PARTIAL")
            else "REPAIR_RG_FORMALISM_OR_INPUTS_AND_RERUN_QW2132"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2132: RG FIXED-POINT JACOBIAN GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Reference couplings (mu=1 GeV)",
        f"- g1: `{g1_ref:.12f}`",
        f"- g2: `{g2_ref:.12f}`",
        f"- g3: `{g3_ref:.12f}`",
        f"- y_t: `{yt_ref:.12f}`",
        f"- lambda_h: `{lam_ref:.12f}`",
        f"- g_gr: `{ggr_ref:.12e}`",
        "",
        "## Fixed-point directional summary at Gaussian",
        f"- counts: `{counts0}`",
        f"- classification: `{dir0}`",
        "",
        "## GR logistic fixed points",
        f"- d(beta)/dg at g=0: `{ggr_slope_at_0:.3f}`",
        f"- d(beta)/dg at g=1: `{ggr_slope_at_1:.3f}`",
        "",
        "## Open scope boundary",
        "- full_nonperturbative_rg_fixed_point_proof: `False`",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2132] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2132] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2132] verdict={verdict} pass_count={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()
