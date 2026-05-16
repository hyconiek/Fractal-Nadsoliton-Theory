#!/usr/bin/env python3
"""P1853 S803 strict B1 symbolic evaluation and scheme trace checkpoint."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def load(name: str) -> dict:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def main() -> None:
    p1850 = load("p1850_s800_strict_gravity_background_b1_symbolic_coefficients_checkpoint.json")
    p1852 = load("p1852_s802_strict_b1_brst_anomaly_and_cutkosky_seed_witness_checkpoint.json")
    p1205 = load("p1205_w4_sympy_cas_runner_summary.json")

    coeffs = ((p1850.get("strict_kernel_to_1loop_symbolic_coefficients") or {}).get("symbolic_coefficients") or {})
    subs = ((p1850.get("strict_kernel_to_1loop_symbolic_coefficients") or {}).get("substitution_tuple") or {})

    sympy_available = bool(p1205.get("sympy_available"))
    values = {}
    eval_trace = {"mode": "symbolic_string_fallback", "sympy_available": sympy_available}

    if sympy_available:
        try:
            import sympy as sp

            omega = sp.nsimplify(subs.get("omega", 0.18575))
            phi = sp.nsimplify(subs.get("phi", 0.16250))
            beta = sp.nsimplify(subs.get("beta", 1.0))
            eta = sp.nsimplify(subs.get("eta", 1.8))
            alpha_geo = 4 * sp.log(2)
            pi = sp.pi

            local = {
                "omega": omega,
                "phi": phi,
                "beta": beta,
                "eta": eta,
                "alpha_geo_strict": alpha_geo,
                "pi": pi,
                "ln": sp.log,
            }
            for k, expr in coeffs.items():
                v = sp.simplify(sp.sympify(expr, locals=local))
                values[k] = {
                    "symbolic": str(v),
                    "numeric_20d": str(sp.N(v, 20)),
                }
            eval_trace = {"mode": "sympy_eval", "sympy_available": True, "sympy_version": sp.__version__}
        except Exception as exc:
            eval_trace = {"mode": "sympy_error_fallback", "sympy_available": True, "error": str(exc)}

    scheme_trace = {
        "scheme_name": "MSbar_B1_seed",
        "epsilon_convention": "d=4-2*epsilon",
        "mu_symbol": "mu_RG",
        "counterterm_identity_ref": "Gamma_1loop_div_B1 + Integral sqrt(-g) * sum_i delta_c_gr_i O_i = 0",
        "normalization_note": "B1 symbolic coefficient layer evaluated in a single declared normalization to avoid cross-scheme drift.",
    }

    out = {
        "packet_id": "P1853",
        "stage_id": "S803",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {
            "p1850_present": "counterterm_cancellation_identity_b1" in p1850,
            "p1852_present": "brst_anomaly_seed_contract" in p1852,
            "p1205_present": "sympy_available" in p1205,
        },
        "b1_symbolic_evaluation": {
            "coefficients_input": coeffs,
            "substitution_tuple": subs,
            "evaluated_coefficients": values,
            "evaluation_trace": eval_trace,
        },
        "renormalization_scheme_trace": scheme_trace,
        "forward_reverse_link": {
            "forward": "K_strict -> symbolic coefficients -> B1 evaluated coefficient layer",
            "reverse": "evaluated coefficient layer -> BRST anomaly seed / Cutkosky seed quantitative backend entry",
            "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        },
        "proven": "B1 symbolic coefficient layer now has a reproducible scheme-tagged backend evaluation trace entry point.",
        "open": "BRST anomaly cancellation coefficients and Cutkosky discontinuity integrals remain uncomputed at theorem level.",
        "false_pass_risk": "Coefficient evaluation trace is not a BRST/Cutkosky/background-independence theorem witness.",
        "next_honest_step": "Compute A_B1 cochain coefficients and first discontinuity integral using the same MSbar_B1_seed normalization and export theorem-linked residual objects.",
        "lay_explanation": "Przepisaliśmy wzory na liczby w jednym ustalonym schemacie, więc łatwiej robić kolejne testy fizyczne, ale to wciąż nie jest końcowy dowód teorii.",
    }

    path = GEN / "p1853_s803_strict_b1_symbolic_evaluation_and_scheme_trace_checkpoint.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
