#!/usr/bin/env python3
"""P1975 S925 strict minimal anisotropic source obligation and energy-sign audit.

P1974 showed a genuine diagonal Bianchi-I metric EOM residual.  This script
constructs the minimal effective source tensor that would cancel that residual
componentwise, then audits its sign.  The result is not a closure proof: the
source is an obligation target unless a strict-side derivation exports it.
"""

from __future__ import annotations

import json
import platform
from pathlib import Path
from typing import Any

import numpy as np
import scipy.linalg as la
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
DEFAULT_TIMESTAMP_UTC = "2026-05-18T00:00:00+00:00"
OUT = GEN / "p1975_s925_strict_minimal_anisotropic_source_obligation_and_energy_sign_audit.json"


def load(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p1974 = load("p1974_s924_strict_bianchi_anisotropic_eom_residual_obstruction_witness.json")

    H, ds1, ds2 = sp.symbols("H dsigma1 dsigma2", real=True)
    s1, s2 = sp.symbols("sigma1 sigma2", real=True)
    local = {"H": H, "sigma1": s1, "sigma2": s2, "dsigma1": ds1, "dsigma2": ds2}

    residual_strings = p1974.get("anisotropic_eom_residual_vector") or []
    residual_vector = [sp.factor(sp.sympify(expr, locals=local)) for expr in residual_strings]
    if len(residual_vector) != 4:
        raise ValueError("P1974 residual vector must have four components")

    q_shear = sp.factor(s1**2 + s1 * s2 + s2**2)
    minimal_delta_t = residual_vector[:]  # In convention G=T, Delta T must equal G_BI-G_FRW.
    post_source_residual = [sp.factor(sp.simplify(r - t)) for r, t in zip(residual_vector, minimal_delta_t)]
    isotropic_subs = {s1: 0, s2: 0, ds1: 0, ds2: 0}
    source_isotropic_limit = [sp.factor(sp.simplify(t.subs(isotropic_subs))) for t in minimal_delta_t]

    rho_req = sp.factor(minimal_delta_t[0])
    spatial = minimal_delta_t[1:]
    mean_pressure_shift = sp.factor(sp.simplify(sum(spatial) / 3))
    traceless_spatial = [sp.factor(sp.simplify(p - mean_pressure_shift)) for p in spatial]
    traceless_sum = sp.factor(sp.simplify(sum(traceless_spatial)))

    # q_shear is positive definite in variables (sigma1, sigma2).  Audit with exact
    # matrix and scipy eigenvalues for machine replay.
    q_matrix = sp.Matrix([[1, sp.Rational(1, 2)], [sp.Rational(1, 2), 1]])
    q_eigs_exact = [sp.factor(ev) for ev in q_matrix.eigenvals().keys()]
    q_eigs_float = la.eigvalsh(np.array(q_matrix.tolist(), dtype=float))
    q_positive_definite = all(float(ev) > 0.0 for ev in q_eigs_float)

    sample_points = [
        {H: sp.Rational(1), s1: sp.Rational(1, 10), s2: sp.Rational(-1, 20), ds1: sp.Rational(1, 100), ds2: sp.Rational(-1, 200)},
        {H: sp.Rational(3, 2), s1: sp.Rational(1, 5), s2: sp.Rational(1, 10), ds1: sp.Rational(-1, 50), ds2: sp.Rational(3, 200)},
        {H: sp.Rational(4, 5), s1: sp.Rational(-1, 8), s2: sp.Rational(1, 16), ds1: sp.Rational(1, 80), ds2: sp.Rational(1, 160)},
    ]
    numeric_rows: list[dict[str, Any]] = []
    for idx, point in enumerate(sample_points):
        delta_vals = [sp.simplify(t.subs(point)) for t in minimal_delta_t]
        post_vals = [sp.simplify(r.subs(point)) for r in post_source_residual]
        rho_val = sp.simplify(rho_req.subs(point))
        q_val = sp.simplify(q_shear.subs(point))
        numeric_rows.append(
            {
                "sample_id": f"minimal_source_sample_{idx}",
                "input": {str(k): str(v) for k, v in point.items()},
                "q_shear": str(q_val),
                "rho_required": str(rho_val),
                "rho_required_float": float(sp.N(rho_val, 30)),
                "minimal_delta_T_components": [str(v) for v in delta_vals],
                "post_source_residual": [str(v) for v in post_vals],
                "post_source_l2_norm": float(la.norm(np.array([float(sp.N(v, 30)) for v in post_vals]), ord=2)),
                "negative_energy_density_for_nonzero_shear": bool(float(sp.N(rho_val, 30)) < 0.0 and float(sp.N(q_val, 30)) > 0.0),
            }
        )

    cancellation_if_admitted = all(expr == 0 for expr in post_source_residual)
    source_frw_limit_zero = all(expr == 0 for expr in source_isotropic_limit)
    rho_equals_minus_q = sp.simplify(rho_req + q_shear) == 0
    sample_cancellation = all(row["post_source_l2_norm"] == 0.0 for row in numeric_rows)
    sample_negative_energy = all(row["negative_energy_density_for_nonzero_shear"] for row in numeric_rows)

    out = {
        "ledger_id": "P1975_S925_STRICT_MINIMAL_ANISOTROPIC_SOURCE_OBLIGATION_AND_ENERGY_SIGN_AUDIT",
        "packet_id": "P1975",
        "stage_id": "S925",
        "produced_by": "p1975_s925_strict_minimal_anisotropic_source_obligation_and_energy_sign_audit.py",
        "timestamp_utc": DEFAULT_TIMESTAMP_UTC,
        "route": "strict_only",
        "depends_on": {
            "p1974_present": p1974.get("result_kind") == "OPEN_OBSTRUCTION_WITH_TRACE",
        },
        "convention": "Einstein component equation G_component = T_component; therefore a cancelling source must satisfy Delta_T = G_BI - G_FRW.",
        "source_obligation": {
            "q_shear": str(q_shear),
            "minimal_delta_T_components": [str(t) for t in minimal_delta_t],
            "rho_required": str(rho_req),
            "mean_spatial_pressure_shift": str(mean_pressure_shift),
            "traceless_spatial_pressure_shift": [str(t) for t in traceless_spatial],
            "traceless_spatial_sum": str(traceless_sum),
            "isotropic_limit": [str(t) for t in source_isotropic_limit],
        },
        "symbolic_checks": {
            "post_source_residual": [str(r) for r in post_source_residual],
            "cancellation_if_source_admitted": cancellation_if_admitted,
            "source_frw_limit_zero": source_frw_limit_zero,
            "rho_required_equals_minus_q_shear": rho_equals_minus_q,
            "q_shear_matrix": str(q_matrix),
            "q_shear_eigenvalues_exact": [str(ev) for ev in q_eigs_exact],
            "q_shear_eigenvalues_scipy": [float(ev) for ev in q_eigs_float],
            "q_shear_positive_definite": q_positive_definite,
        },
        "numeric_probe_table": numeric_rows,
        "gatekeeper_checks": {
            "minimal_source_cancels_residual_if_admitted": cancellation_if_admitted and sample_cancellation,
            "minimal_source_has_zero_frw_limit": source_frw_limit_zero,
            "required_energy_density_negative_for_nonzero_shear": rho_equals_minus_q and q_positive_definite and sample_negative_energy,
            "strict_source_derivation_exported": False,
        },
        "result_kind": "OPEN_OBSTRUCTION_WITH_TRACE",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "obstruction_tags": [
            "MINIMAL_SOURCE_IS_OBLIGATION_NOT_DERIVED_STRICT_SOURCE",
            "REQUIRED_RHO_EQUALS_NEGATIVE_SHEAR_QUADRATIC",
            "POSITIVE_ENERGY_STRICT_SOURCE_NOT_AVAILABLE",
            "BACKGROUND_INDEPENDENCE_STILL_OPEN",
        ],
        "theorem_export": {
            "conditional_positive_statement": "If a strict source with Delta_T equal to the exported minimal components is admitted, the P1974 FRW/Bianchi-I component residual cancels exactly and the source vanishes in the FRW isotropic limit.",
            "negative_energy_statement": "Under the convention G=T, the required energy component is rho_required=-(sigma1^2+sigma1*sigma2+sigma2^2), which is strictly negative for nonzero trace-free shear because the quadratic form has positive eigenvalues 1/2 and 3/2.",
            "missing_strict_provider": "No strict-side derivation of this source is exported in the current repo state; admitting it as an axiom would be non-strict.",
            "not_licensed": [
                "positive-energy anisotropic source closure",
                "global background-independence closure",
                "PO2/PO3 closure",
                "BRST or Cutkosky closure",
                "QW-2191 selector discharge",
                "ToE closure",
            ],
        },
        "false_pass_guard": "P1975 constructs the exact source obligation and a sign obstruction; it does not derive the source from K_strict or L_total and therefore cannot close background-independence.",
        "next_honest_step": "Audit the full strict L_total sector list for a derived shear-energy term with the required negative rho signature, or prove a no-go theorem that no positive-energy strict anisotropic provider can cancel P1974.",
        "lay_explanation": "Wyliczyliśmy dokładnie, jakie dodatkowe źródło musiałoby istnieć, żeby usunąć anizotropowy błąd z P1974.  Matematycznie takie źródło skasowałoby błąd, ale jego energia wychodzi ujemna dla niezerowej anizotropii, więc nie wolno go po prostu dopisać jako fizycznego składnika teorii bez osobnego ścisłego wyprowadzenia.",
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": __import__("scipy").__version__,
            "sympy": sp.__version__,
        },
    }
    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
