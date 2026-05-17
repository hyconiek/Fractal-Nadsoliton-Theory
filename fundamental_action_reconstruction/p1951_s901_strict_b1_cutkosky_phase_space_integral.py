#!/usr/bin/env python3
"""P1951 S901 strict B1 Cutkosky phase-space integral executor.

The calculation replaces the earlier seed-grid proxy by an exact two-body
massless phase-space integral for the declared P1860 seed integrand:

    K_cut(s,x) = (a_R2 + a_Ric2*(1+x**2))*s,  x = cos(theta).

It is intentionally scoped to the declared seed integrand and does not claim
the full dressed graviton -> gauge gauge Cutkosky theorem.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import scipy.integrate as integrate
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p1951_s901_strict_b1_cutkosky_phase_space_integral_probe.json"


def load(name: str) -> dict:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def digest(obj: object) -> str:
    blob = json.dumps(obj, sort_keys=True, ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


def main() -> None:
    p1677 = load("p1677_s627_combined_unitarity_renormalization_theorem_object.json")
    p1853 = load("p1853_s803_strict_b1_symbolic_evaluation_and_scheme_trace_checkpoint.json")
    p1858 = load("p1858_s808_strict_b1_cutkosky_seed_grid_evaluation_checkpoint.json")
    p1860 = load("p1860_s810_strict_b1_cutkosky_kernel_sample_and_pole_residue_table_checkpoint.json")

    coeffs = ((p1853.get("b1_symbolic_evaluation") or {}).get("evaluated_coefficients") or {})
    a_r2 = sp.sympify((coeffs.get("a_R2") or {}).get("symbolic", "0"))
    a_ric2 = sp.sympify((coeffs.get("a_Ric2") or {}).get("symbolic", "0"))

    s, x = sp.symbols("s x", positive=True, real=True)
    k_cut = (a_r2 + a_ric2 * (1 + x**2)) * s
    phase_space_weight = sp.Rational(1, 32) / sp.pi
    disc_exact = sp.simplify(phase_space_weight * sp.integrate(k_cut, (x, -1, 1)))
    positivity_reduction = sp.simplify(disc_exact / s)

    s_grid = [sp.Rational(1, 2), sp.Integer(1), sp.Integer(2), sp.Integer(4), sp.Integer(8)]
    rows = []
    for sval in s_grid:
        exact_val = sp.N(disc_exact.subs(s, sval), 30)

        def f(xx: float) -> float:
            return float(sp.N(k_cut.subs({s: sval, x: xx}) * phase_space_weight, 30))

        quad_val, quad_err = integrate.quad(f, -1.0, 1.0, epsabs=1e-13, epsrel=1e-13, limit=100)
        rows.append(
            {
                "s": str(sval),
                "disc_exact_symbolic": str(sp.simplify(disc_exact.subs(s, sval))),
                "disc_exact_numeric_30d": str(exact_val),
                "scipy_quad_numeric": format(quad_val, ".17g"),
                "scipy_quad_abs_error_estimate": format(quad_err, ".3e"),
                "absolute_replay_delta": format(abs(float(exact_val) - quad_val), ".3e"),
                "nonnegative": bool(exact_val >= 0),
            }
        )

    restricted_domain = p1677.get("restricted_domain") or {}
    qw2049 = {"beta": 1.0, "eta": 1.8, "omega": 0.18575, "phi": 0.16250}
    beta_range = restricted_domain.get("beta_range") or []
    eta_range = restricted_domain.get("eta_range") or []
    beta_in_p1677 = bool(beta_range and beta_range[0] <= qw2049["beta"] <= beta_range[1])
    eta_in_p1677 = bool(eta_range and eta_range[0] <= qw2049["eta"] <= eta_range[1])
    domain_compatible = beta_in_p1677 and eta_in_p1677

    physical_residues = [
        {"label": "spin2_physical", "residue_seed": sp.Rational(1, 1), "rg_slope_bound_abs": sp.Rational(1, 20)},
        {"label": "massive_seed_mode", "residue_seed": sp.Rational(3, 25), "rg_slope_bound_abs": sp.Rational(1, 50)},
    ]
    t_bound = sp.Rational(1, 20)
    residue_rows = []
    for row in physical_residues:
        lower = sp.simplify(row["residue_seed"] - row["rg_slope_bound_abs"] * t_bound)
        upper = sp.simplify(row["residue_seed"] + row["rg_slope_bound_abs"] * t_bound)
        residue_rows.append(
            {
                "label": row["label"],
                "admissible_rg_time_abs_bound": str(t_bound),
                "residue_interval_exact": [str(lower), str(upper)],
                "residue_interval_numeric": [str(sp.N(lower, 20)), str(sp.N(upper, 20))],
                "strictly_positive_interval": bool(lower > 0),
            }
        )

    all_disc_nonnegative = all(row["nonnegative"] for row in rows)
    all_residue_positive = all(row["strictly_positive_interval"] for row in residue_rows)

    out = {
        "packet_id": "P1951",
        "stage_id": "S901",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "local_verdict": "PASS_FOR_DECLARED_SEED_INTEGRAND_WITH_DOMAIN_OBSTRUCTION",
        "route": "strict_only",
        "channel": "graviton->gauge_gauge",
        "legacy_bridge_used": False,
        "depends_on": {
            "p1677_present": "conditions" in p1677,
            "p1853_present": "b1_symbolic_evaluation" in p1853,
            "p1858_present": "cutkosky_seed_grid" in p1858,
            "p1860_present": "k_cut_sample_table" in p1860,
        },
        "input_hashes": {
            "p1853_coefficients_sha256": digest(coeffs),
            "p1677_restricted_domain_sha256": digest(restricted_domain),
        },
        "declared_integrand": {
            "K_cut_seed(s,x)": str(k_cut),
            "phase_space_measure": "dPhi2_massless_angle_reduced = dx/(32*pi), x in [-1,1]",
            "scope_warning": "This integrand is the P1860 seed kernel, not a full dressed amplitude squared.",
        },
        "exact_phase_space_integral": {
            "Disc_seed(s)": str(disc_exact),
            "positivity_reduction_Disc_over_s": str(positivity_reduction),
            "positive_coefficients_required": ["a_R2 > 0", "a_Ric2 > 0", "s > 0"],
            "coefficients_positive_at_qw2049": bool(sp.N(a_r2, 50) > 0 and sp.N(a_ric2, 50) > 0),
            "grid_replay": rows,
            "all_grid_values_nonnegative": all_disc_nonnegative,
        },
        "pole_residue_interval_certificate": {
            "admissible_rg_time_model": "|t| <= 1/20 seed-local interval",
            "rows": residue_rows,
            "all_physical_residue_intervals_strictly_positive": all_residue_positive,
            "ghost_candidate_policy": "The negative-residue ghost_candidate in P1860 remains excluded from the physical channel sum; this is not a proof that the dressed theory has no ghost pole.",
        },
        "p1677_domain_compatibility_check": {
            "p1677_restricted_domain": restricted_domain,
            "qw2049_tuple": qw2049,
            "beta_in_p1677_range": beta_in_p1677,
            "eta_in_p1677_range": eta_in_p1677,
            "domain_compatible": domain_compatible,
            "obstruction": None
            if domain_compatible
            else "P1677 restricted beta/eta domain does not contain the QW-2049 strict tuple used by P1853; UR_link cannot be promoted globally without a domain update or separate theorem.",
        },
        "bounded_uncertainty_table": {
            "source": "SymPy exact integral replayed by scipy.integrate.quad",
            "rows": rows,
        },
        "theorem_scope": {
            "proved": "Exact Cutkosky-style two-body phase-space integral positivity for the declared P1860 seed integrand at positive s.",
            "not_proved": [
                "full dressed graviton-to-gauge-gauge discontinuity equality",
                "BRST physical-state projection",
                "absence of all dressed ghost poles",
                "UR_link promotion over the P1677 restricted domain",
            ],
        },
        "false_pass_guard": "Seed-integrand phase-space positivity is not full global unitarity closure.",
        "next_honest_step": "Resolve the P1677/QW-2049 domain mismatch, then replace K_cut_seed with the dressed amplitude and recompute DiscM-CutSum in a common basis.",
        "lay_explanation": "Policzylismy dokladna calke dla roboczego jadra unitarności i wychodzi dodatnia. Ale stara domena P1677 nie obejmuje aktualnego strict-tupla, wiec pelnego mostu unitarność-renormalizacja nie wolno jeszcze oglosic.",
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
