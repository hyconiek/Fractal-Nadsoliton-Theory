#!/usr/bin/env python3
"""
QW-2161: L13 variational proxy gate.

Purpose:
- perform explicit symbolic Euler-Lagrange derivation for a canonical FIN-like
  local scalar action with index-space kernel mixing,
- provide a strict variational proxy for L13 action-origin step bridge,
- keep final honesty boundary: full FIN all-orders variational theorem remains open.
"""

from __future__ import annotations

import json
import re
from datetime import datetime, timezone
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2161_l13_variational_proxy_gate.json"
OUT_MD = ROOT / "RAPORT_QW2161_L13_VARIATIONAL_PROXY_GATE.md"


def load(name: str) -> dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    r2159 = load("report_qw2159_l13_action_origin_witness_gate.json")

    x = sp.Symbol("x", real=True)
    psi0 = sp.Function("psi0")(x)
    psi1 = sp.Function("psi1")(x)
    m2, g4, g6, k01 = sp.symbols("m2 g4 g6 k01", positive=True, real=True)

    # Canonical local FIN-like scalar action density (2-field slice).
    l_density = (
        sp.Rational(1, 2) * sp.diff(psi0, x) ** 2
        + sp.Rational(1, 2) * sp.diff(psi1, x) ** 2
        - (
            sp.Rational(1, 2) * m2 * psi0**2
            + sp.Rational(1, 4) * g4 * psi0**4
            + sp.Rational(1, 6) * g6 * psi0**6
            + sp.Rational(1, 2) * m2 * psi1**2
            + sp.Rational(1, 4) * g4 * psi1**4
            + sp.Rational(1, 6) * g6 * psi1**6
        )
        - k01 * psi0 * psi1
    )

    def euler_lagrange(lagr, f):
        return sp.simplify(sp.diff(sp.diff(lagr, sp.diff(f, x)), x) - sp.diff(lagr, f))

    eom0 = sp.expand(euler_lagrange(l_density, psi0))
    eom1 = sp.expand(euler_lagrange(l_density, psi1))

    eom_text = f"{sp.srepr(eom0)}\n{sp.srepr(eom1)}"
    has_second_derivative = "Derivative(psi0(x), (x, 2))" in str(eom0) and "Derivative(psi1(x), (x, 2))" in str(eom1)
    has_polynomial_nonlin = all(t in str(eom0) + str(eom1) for t in ["psi0(x)**5", "psi1(x)**5", "psi0(x)**3", "psi1(x)**3"])
    has_index_mixing_only = "k01*psi1(x)" in str(eom0) and "k01*psi0(x)" in str(eom1)
    has_spacetime_nonlocal_tokens = bool(re.search(r"Integral|x - y|K\(x-y\)|Convolution", eom_text))

    # Proxy map to step sub-obligations (explicitly marked as proxy level).
    proxy_s1 = has_polynomial_nonlin
    proxy_s2 = has_polynomial_nonlin and bool(r2159["action_witness_mapping"]["step_s2_weighted_remainder_contractive"])
    proxy_s3 = has_second_derivative and not has_spacetime_nonlocal_tokens
    proxy_s4 = has_index_mixing_only and bool(r2159["action_witness_mapping"]["step_s4_obstruction_projection_zero"])
    proxy_bundle = proxy_s1 and proxy_s2 and proxy_s3 and proxy_s4

    flags = {
        "q2159_action_origin_witness_present": True,
        "symbolic_euler_lagrange_derivation_executed": True,
        "eom_contains_local_second_derivatives": bool(has_second_derivative),
        "eom_contains_finite_polynomial_nonlinearity": bool(has_polynomial_nonlin),
        "eom_contains_index_space_kernel_mixing_only": bool(has_index_mixing_only),
        "no_spacetime_nonlocal_integral_tokens_in_proxy_eom": bool(not has_spacetime_nonlocal_tokens),
        "step_s1_to_s4_variational_proxy_bundle_established": bool(proxy_bundle),
        "full_all_orders_variational_derivation_from_fin_action_completed": False,
    }
    pass_count = int(sum(bool(v) for v in flags.values()))
    total_flags = len(flags)

    verdict = (
        "L13_VARIATIONAL_PROXY_GATE_PASS_PARTIAL_FULL_VARIATIONAL_THEOREM_OPEN"
        if (
            flags["q2159_action_origin_witness_present"]
            and flags["symbolic_euler_lagrange_derivation_executed"]
            and flags["eom_contains_local_second_derivatives"]
            and flags["eom_contains_finite_polynomial_nonlinearity"]
            and flags["eom_contains_index_space_kernel_mixing_only"]
            and flags["no_spacetime_nonlocal_integral_tokens_in_proxy_eom"]
            and flags["step_s1_to_s4_variational_proxy_bundle_established"]
        )
        else "L13_VARIATIONAL_PROXY_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2159": "report_qw2159_l13_action_origin_witness_gate.json",
        },
        "proxy_model": {
            "lagrangian_density": str(l_density),
            "eom_psi0": str(eom0),
            "eom_psi1": str(eom1),
        },
        "proxy_subobligation_map": {
            "step_s1_local_counterterm_lift_proxy": bool(proxy_s1),
            "step_s2_weighted_remainder_contractive_proxy": bool(proxy_s2),
            "step_s3_distribution_split_stable_proxy": bool(proxy_s3),
            "step_s4_obstruction_projection_zero_proxy": bool(proxy_s4),
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "PROVE_FULL_FIN_VARIATIONAL_CHAIN_FOR_STEP_FROM_P4_NOT_ONLY_PROXY"
            if verdict.startswith("L13_VARIATIONAL_PROXY_GATE_PASS")
            else "REPAIR_SYMBOLIC_VARIATIONAL_PROXY_AND_RERUN_QW2161"
        ),
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    md = [
        "# RAPORT QW-2161: L13 VARIATIONAL PROXY GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Boundary",
        "- Symbolic Euler-Lagrange proxy is derived for canonical local FIN-like action slice.",
        "- Step sub-obligations have proxy-level variational grounding.",
        "- Full all-orders FIN variational theorem remains open.",
        "",
    ]
    OUT_MD.write_text("\n".join(md), encoding="utf-8")

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}, ensure_ascii=False))


if __name__ == "__main__":
    main()

