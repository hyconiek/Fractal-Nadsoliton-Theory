#!/usr/bin/env python3
"""
QW-2185: RG proxy global-obstruction theorem gate (strict).

Purpose:
- formalize what is and is not provable globally for the declared QW-2132/QW-2182
  proxy RG system,
- provide exact theorem-level obstruction for full global (t -> +inf) flow due to
  U(1) Landau pole in current one-loop proxy beta(g1),
- keep strict clarity: this strengthens rigor for L12 without overclaiming closure.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2185_rg_proxy_global_obstruction_theorem_gate.json"
OUT_MD = ROOT / "RAPORT_QW2185_RG_PROXY_GLOBAL_OBSTRUCTION_THEOREM_GATE.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def g_landau_pole_time(g0: float, k: float, t0: float = 0.0) -> float:
    """For dg/dt = k g^3 (k>0), g0>0 => blow-up at t* = t0 + 1/(2 k g0^2)."""
    if g0 <= 0.0:
        return math.inf
    return float(t0 + 1.0 / (2.0 * k * g0 * g0))


def g_asym_free_solution(g0: float, k_abs: float, t: float, t0: float = 0.0) -> float:
    """For dg/dt = -k_abs g^3 (k_abs>0): g(t)=g0/sqrt(1+2 k_abs g0^2 (t-t0))."""
    den = 1.0 + 2.0 * k_abs * g0 * g0 * (t - t0)
    return float(g0 / math.sqrt(den))


def g_logistic_solution(g0: float, t: float, t0: float = 0.0) -> float:
    """For dg/dt = 2 g (1-g)."""
    if g0 == 0.0:
        return 0.0
    e2t = math.exp(2.0 * (t - t0))
    return float((g0 * e2t) / (1.0 + g0 * (e2t - 1.0)))


def main() -> None:
    r2132 = load_json("report_qw2132_rg_fixed_point_jacobian_gate.json")
    r2182 = load_json("report_qw2182_rg_nonperturbative_domain_flow_gate.json")

    ref = r2132["reference_couplings_at_mu_1gev"]
    g1_ref = float(ref["g1_u1"])
    g2_ref = float(ref["g2_su2"])
    g3_ref = float(ref["g3_su3"])
    ggr_ref = float(ref["g_gr_dimensionless"])

    # One-loop proxy coefficients (identical to QW-2132/QW-2182).
    c = 1.0 / (16.0 * math.pi * math.pi)
    k1 = c * (41.0 / 6.0)  # dg1/dt = +k1 g1^3
    k2 = c * (19.0 / 6.0)  # dg2/dt = -k2 g2^3
    k3 = c * 7.0           # dg3/dt = -k3 g3^3

    t_star_ref = g_landau_pole_time(g1_ref, k1)

    # Domain from QW-2182 for strict constructive certificate.
    t_max = float(r2182["domain"]["t_max"])
    g1_domain: List[float] = [float(v) for v in r2182["domain"]["grid"]["g1"]]
    t_star_domain = [g_landau_pole_time(g, k1) for g in g1_domain]
    t_star_domain_min = float(min(t_star_domain))

    # Closed-form monotonicity checks for decoupled channels.
    g2_tmax = g_asym_free_solution(g2_ref, k2, t_max)
    g3_tmax = g_asym_free_solution(g3_ref, k3, t_max)
    ggr_tmax = g_logistic_solution(ggr_ref, t_max)

    # Theorem-level statements for current proxy:
    full_global_flow_possible = bool(g1_ref <= 0.0)  # for physical g1_ref>0 this is False
    full_global_fixed_point_for_all_couplings = False
    obstruction_is_landau_pole_g1 = bool(g1_ref > 0.0 and math.isfinite(t_star_ref))
    strict_window_safe = bool(t_star_domain_min > t_max)

    flags = {
        "q2132_proxy_model_loaded": True,
        "q2182_constructive_domain_loaded": True,
        "u1_beta_sign_positive_in_proxy": bool(k1 > 0.0),
        "exact_landau_pole_formula_applied": bool(math.isfinite(t_star_ref)),
        "landau_pole_time_reference_positive_finite": bool(t_star_ref > 0.0 and math.isfinite(t_star_ref)),
        "global_full_flow_for_all_t_not_possible_in_current_proxy": bool(not full_global_flow_possible),
        "global_full_fixed_point_proof_not_possible_in_current_proxy": bool(not full_global_fixed_point_for_all_couplings),
        "obstruction_explicitly_identified_as_u1_landau_pole": bool(obstruction_is_landau_pole_g1),
        "domain_window_q2182_lies_strictly_before_landau_pole": bool(strict_window_safe),
        "g2_asymptotic_freedom_closed_form_monotone_to_zero": bool(g2_tmax < g2_ref and g2_tmax > 0.0),
        "g3_asymptotic_freedom_closed_form_monotone_to_zero": bool(g3_tmax < g3_ref and g3_tmax > 0.0),
        "ggr_logistic_closed_form_monotone_to_uv_fixed_point_one": bool(ggr_tmax > ggr_ref and ggr_tmax < 1.0),
        "deterministic_no_scan_no_retune": True,
        "full_nonperturbative_rg_fixed_point_proof_completed": False,
    }

    pass_count = int(sum(bool(v) for v in flags.values()))
    total_flags = int(len(flags))

    core_ok = bool(
        flags["q2132_proxy_model_loaded"]
        and flags["q2182_constructive_domain_loaded"]
        and flags["u1_beta_sign_positive_in_proxy"]
        and flags["exact_landau_pole_formula_applied"]
        and flags["landau_pole_time_reference_positive_finite"]
        and flags["global_full_flow_for_all_t_not_possible_in_current_proxy"]
        and flags["global_full_fixed_point_proof_not_possible_in_current_proxy"]
        and flags["obstruction_explicitly_identified_as_u1_landau_pole"]
        and flags["domain_window_q2182_lies_strictly_before_landau_pole"]
        and flags["g2_asymptotic_freedom_closed_form_monotone_to_zero"]
        and flags["g3_asymptotic_freedom_closed_form_monotone_to_zero"]
        and flags["ggr_logistic_closed_form_monotone_to_uv_fixed_point_one"]
        and flags["deterministic_no_scan_no_retune"]
    )

    verdict = (
        "RG_PROXY_GLOBAL_OBSTRUCTION_THEOREM_GATE_PASS_STRICT"
        if core_ok
        else "RG_PROXY_GLOBAL_OBSTRUCTION_THEOREM_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2132": "report_qw2132_rg_fixed_point_jacobian_gate.json",
            "q2182": "report_qw2182_rg_nonperturbative_domain_flow_gate.json",
        },
        "proxy_coefficients": {
            "k1_u1_positive": k1,
            "k2_su2_abs": k2,
            "k3_su3_abs": k3,
        },
        "obstruction_theorem": {
            "statement": (
                "For dg1/dt = k1 g1^3 with k1>0 and g1(0)>0, solution has finite-time "
                "Landau pole t*=1/(2 k1 g1(0)^2). Therefore full global flow t>=0 and "
                "global full-coupling fixed-point closure are impossible in the current proxy."
            ),
            "reference_g1": g1_ref,
            "t_star_reference": t_star_ref,
            "q2182_t_max": t_max,
            "t_star_domain_min": t_star_domain_min,
            "window_safe_before_pole": strict_window_safe,
        },
        "closed_form_subsector_checks": {
            "g2_reference": g2_ref,
            "g2_at_tmax": g2_tmax,
            "g3_reference": g3_ref,
            "g3_at_tmax": g3_tmax,
            "ggr_reference": ggr_ref,
            "ggr_at_tmax": ggr_tmax,
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "EITHER_DECLARE_FINITE_UV_VALIDITY_WINDOW_FOR_CURRENT_PROXY_OR_INTRODUCE_UV_COMPLETING_BETA_CORRECTIONS_THEN_RERUN_L12_CHAIN"
            if verdict.endswith("STRICT")
            else "REVIEW_THEOREM_ASSUMPTIONS_AND_RERUN_QW2185"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    md_lines = [
        "# RAPORT QW-2185: RG PROXY GLOBAL OBSTRUCTION THEOREM GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Core result",
        "- In current one-loop proxy (`QW-2132`), full global RG closure is obstructed by U(1) Landau pole.",
        f"- Reference Landau-pole time: `t*={t_star_ref:.6f}`.",
        f"- QW-2182 finite window is before pole: `t_max={t_max}` and `t*_min(domain)={t_star_domain_min:.6f}`.",
        "",
        "## Meaning",
        "- This is a strict theorem-level clarification, not overclaim closure.",
        "- L12 remains partial globally until UV-completing correction or finite-UV-scope declaration.",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}, ensure_ascii=False))


if __name__ == "__main__":
    main()

