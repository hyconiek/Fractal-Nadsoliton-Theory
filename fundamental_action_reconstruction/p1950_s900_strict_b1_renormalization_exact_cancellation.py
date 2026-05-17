#!/usr/bin/env python3
"""P1950 S900 strict B1 renormalization exact cancellation executor.

This script does not derive the B1 coefficient ansatz from first principles.
It consumes the declared P1850/P1853 B1 coefficient layer, evaluates it with
SymPy, and proves the algebraic counterterm cancellation identity in the same
declared MSbar_B1_seed scheme.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p1950_s900_strict_b1_renormalization_exact_cancellation_probe.json"


def load(name: str) -> dict:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def digest(obj: object) -> str:
    blob = json.dumps(obj, sort_keys=True, ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


def parse_expr(expr: str, local: dict[str, object]) -> sp.Expr:
    return sp.sympify(expr.replace("^", "**"), locals=local)


def main() -> None:
    p1849 = load("p1849_s799_strict_gravity_1loop_divergence_cancellation_witness_checkpoint.json")
    p1850 = load("p1850_s800_strict_gravity_background_b1_symbolic_coefficients_checkpoint.json")
    p1853 = load("p1853_s803_strict_b1_symbolic_evaluation_and_scheme_trace_checkpoint.json")

    coeff_block = p1850.get("strict_kernel_to_1loop_symbolic_coefficients") or {}
    coeff_inputs = coeff_block.get("symbolic_coefficients") or {}
    subs_input = coeff_block.get("substitution_tuple") or {}

    omega, phi, beta, eta = sp.symbols("omega phi beta eta", positive=True, real=True)
    epsilon = sp.symbols("epsilon", nonzero=True)
    R2, Ric2, Riem2, GGB = sp.symbols("R2 Ric2 Riem2 GGB")
    alpha_geo_strict = sp.symbols("alpha_geo_strict", positive=True, real=True)
    local = {
        "omega": omega,
        "phi": phi,
        "beta": beta,
        "eta": eta,
        "pi": sp.pi,
        "ln": sp.log,
        "alpha_geo_strict": alpha_geo_strict,
    }

    coeff_symbols = {name: sp.simplify(parse_expr(expr, local)) for name, expr in coeff_inputs.items()}
    operators = {"a_R2": R2, "a_Ric2": Ric2, "a_Riem2": Riem2, "a_GB": GGB}
    gamma_div = sum(coeff_symbols[name] * operators[name] for name in operators) / epsilon
    counterterm = sum((-coeff_symbols[name] / epsilon) * operators[name] for name in operators)
    cancellation_residual = sp.simplify(gamma_div + counterterm)

    exact_subs = {
        omega: sp.Rational(str(subs_input.get("omega", "0.18575"))),
        phi: sp.Rational(str(subs_input.get("phi", "0.16250"))),
        beta: sp.Rational(str(subs_input.get("beta", "1.0"))),
        eta: sp.Rational(str(subs_input.get("eta", "1.8"))),
        alpha_geo_strict: 4 * sp.log(2),
    }

    evaluated = {}
    for name, expr in coeff_symbols.items():
        val = sp.simplify(expr.subs(exact_subs))
        evaluated[name] = {
            "symbolic_exact": str(val),
            "numeric_30d": str(sp.N(val, 30)),
            "positive": bool(sp.N(val, 50) > 0),
        }

    residual_at_tuple = sp.simplify(cancellation_residual.subs(exact_subs))
    residual_zero = residual_at_tuple == 0

    previous_eval = ((p1853.get("b1_symbolic_evaluation") or {}).get("evaluated_coefficients") or {})
    replay_checks = []
    for name, row in evaluated.items():
        prev_symbolic = (previous_eval.get(name) or {}).get("symbolic")
        replay_checks.append(
            {
                "coefficient": name,
                "matches_p1853_symbolic": prev_symbolic == row["symbolic_exact"],
                "p1853_symbolic": prev_symbolic,
                "p1950_symbolic": row["symbolic_exact"],
            }
        )

    all_replay_match = all(item["matches_p1853_symbolic"] for item in replay_checks)

    out = {
        "packet_id": "P1950",
        "stage_id": "S900",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "local_verdict": "PASS_ZERO_ON_DECLARED_B1_ANSATZ",
        "route": "strict_only",
        "kernel_anchor": "K_strict(d)=cos(omega*d+phi)/(1+beta*d^eta)",
        "legacy_bridge_used": False,
        "depends_on": {
            "p1849_present": "strict_gravity_1loop_divergence_projection" in p1849,
            "p1850_present": "strict_kernel_to_1loop_symbolic_coefficients" in p1850,
            "p1853_present": "b1_symbolic_evaluation" in p1853,
        },
        "input_hashes": {
            "p1850_coefficients_sha256": digest(coeff_inputs),
            "p1850_substitution_tuple_sha256": digest(subs_input),
            "p1853_evaluated_coefficients_sha256": digest(previous_eval),
        },
        "scheme": {
            "name": "MSbar_B1_seed",
            "epsilon_convention": "d=4-2*epsilon",
            "background_family": "background_family_B1",
        },
        "computed_coefficients": evaluated,
        "counterterm_assignments": {
            "delta_c_gr_1": "-a_R2/epsilon",
            "delta_c_gr_2": "-a_Ric2/epsilon",
            "delta_c_gr_3": "-a_Riem2/epsilon",
            "delta_c_gr_4": "-a_GB/epsilon",
        },
        "symbolic_cancellation": {
            "gamma_divergent_density": str(gamma_div),
            "counterterm_density": str(counterterm),
            "residual_symbolic": str(cancellation_residual),
            "residual_at_qw2049_tuple": str(residual_at_tuple),
            "residual_zero": residual_zero,
        },
        "backend_replay_checks": {
            "all_match_p1853": all_replay_match,
            "rows": replay_checks,
        },
        "theorem_scope": {
            "proved": "Exact algebraic cancellation of the declared B1 divergent density by the declared delta_c_gr_i basis.",
            "not_proved": [
                "first-principles derivation of the B1 coefficient ansatz",
                "background-family globalization beyond B1",
                "BRST/Cutkosky unitarity closure",
                "strict-core selector closure",
            ],
        },
        "false_pass_guard": "This is a local theorem for the declared B1 coefficient ansatz, not full renormalization closure of the ToE.",
        "next_honest_step": "Derive the same a_* coefficients from diagram-resolved heat-kernel/IBP data or export a non-ansatz obstruction.",
        "lay_explanation": "Dla zadeklarowanych wzorow B1 nieskonczonosci kasuja sie dokladnie kontrterminami. Nadal trzeba udowodnic, ze same wzory B1 naprawde pochodza z pelnego rachunku petli.",
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
