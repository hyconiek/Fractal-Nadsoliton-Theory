#!/usr/bin/env python3
"""P1978 S928 strict energy-neutral tensor-transport obstruction.

This packet takes the next honest escape route after P1977 and tests whether a
non-minimal tensorial transport connection can remove the P1974 Bianchi-I
residual without introducing an energy-density component.  The result is a
bounded negative theorem: any energy-neutral spatial/shear-only transport leaves
the G00 residual -Q_shear untouched, and any trace-free spatial transport also
cannot cancel the spatial trace -3 Q_shear.

Scope: this does not rule out a future strict gravitational/shear sector with a
derived lapse/energy contribution, nor does it close background independence.
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
OUT = GEN / "p1978_s928_strict_energy_neutral_tensor_transport_obstruction.json"


def load(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p1974 = load("p1974_s924_strict_bianchi_anisotropic_eom_residual_obstruction_witness.json")
    p1977 = load("p1977_s927_strict_positive_energy_anisotropic_provider_bounded_no_go.json")

    H, sigma1, sigma2, dsigma1, dsigma2 = sp.symbols("H sigma1 sigma2 dsigma1 dsigma2", real=True)
    sigma3 = -sigma1 - sigma2
    dsigma3 = -dsigma1 - dsigma2
    q_shear = sp.factor(sigma1**2 + sigma1 * sigma2 + sigma2**2)
    s1 = sp.factor(3 * H * sigma1 + dsigma1)
    s2 = sp.factor(3 * H * sigma2 + dsigma2)
    s3 = sp.factor(3 * H * sigma3 + dsigma3)
    residual = [sp.factor(-q_shear), sp.factor(s1 - q_shear), sp.factor(s2 - q_shear), sp.factor(s3 - q_shear)]

    # Energy-neutral non-minimal transport ansatz.  The spatial slots u_i are
    # intentionally unconstrained first; trace-free is checked as a stricter
    # shear-only subcase.
    u0, u1, u2, u3 = sp.symbols("u0 u1 u2 u3", real=True)
    transport = [u0, u1, u2, u3]
    post_transport = [sp.factor(r - t) for r, t in zip(residual, transport)]

    energy_neutral_condition = sp.Eq(u0, 0)
    required_u0_for_full_cancellation = sp.solve(sp.Eq(post_transport[0], 0), u0)[0]
    energy_neutral_energy_residual = sp.factor(post_transport[0].subs(u0, 0))
    # Compatibility of U00=0 with required_U00=-Q_shear is the zero-set
    # condition Q_shear=0.  Keep it explicit; no branch exemplar is used.
    energy_neutral_requires_q_zero = sp.factor(-required_u0_for_full_cancellation)

    spatial_trace_residual = sp.factor(sum(residual[1:]))
    spatial_trace_transport = sp.factor(u1 + u2 + u3)
    post_spatial_trace = sp.factor(spatial_trace_residual - spatial_trace_transport)
    tracefree_post_spatial_trace = sp.factor(post_spatial_trace.subs(u3, -u1 - u2))

    q_matrix = sp.Matrix([[1, sp.Rational(1, 2)], [sp.Rational(1, 2), 1]])
    scipy_eigs = la.eigvalsh(np.array(q_matrix.tolist(), dtype=float))
    q_positive_definite = all(float(ev) > 0 for ev in scipy_eigs)

    sample_points = [
        {H: sp.Rational(1, 1), sigma1: sp.Rational(1, 10), sigma2: sp.Rational(-1, 20), dsigma1: sp.Rational(1, 100), dsigma2: sp.Rational(-1, 200)},
        {H: sp.Rational(3, 2), sigma1: sp.Rational(1, 5), sigma2: sp.Rational(1, 10), dsigma1: sp.Rational(-1, 50), dsigma2: sp.Rational(3, 200)},
        {H: sp.Rational(4, 5), sigma1: sp.Rational(-1, 8), sigma2: sp.Rational(1, 16), dsigma1: sp.Rational(1, 80), dsigma2: sp.Rational(1, 160)},
        {H: sp.Rational(5, 4), sigma1: sp.Rational(2, 7), sigma2: sp.Rational(-1, 9), dsigma1: sp.Rational(1, 70), dsigma2: sp.Rational(-1, 90)},
    ]
    rows: list[dict[str, Any]] = []
    for idx, point in enumerate(sample_points):
        q_val = sp.simplify(q_shear.subs(point))
        residual_val = [sp.simplify(expr.subs(point)) for expr in residual]
        energy_neutral_leftover = sp.simplify(energy_neutral_energy_residual.subs(point))
        tracefree_spatial_leftover = sp.simplify(tracefree_post_spatial_trace.subs(point))
        rows.append(
            {
                "sample_id": f"energy_neutral_transport_sample_{idx}",
                "input": {str(k): str(v) for k, v in point.items()},
                "q_shear": str(q_val),
                "residual_vector": [str(v) for v in residual_val],
                "energy_neutral_leftover_G00": str(energy_neutral_leftover),
                "tracefree_spatial_leftover_sum": str(tracefree_spatial_leftover),
                "energy_neutral_full_cancellation_impossible": bool(float(sp.N(q_val, 30)) > 0.0 and float(sp.N(energy_neutral_leftover, 30)) < 0.0),
                "tracefree_spatial_cancellation_impossible": bool(float(sp.N(q_val, 30)) > 0.0 and float(sp.N(tracefree_spatial_leftover, 30)) < 0.0),
            }
        )

    p1974_present = p1974.get("result_kind") == "OPEN_OBSTRUCTION_WITH_TRACE"
    p1977_no_go = p1977.get("result_kind") == "PASS_BOUNDED_NO_GO"
    energy_symbolic_obstruction = bool(energy_neutral_energy_residual == -q_shear and q_positive_definite)
    tracefree_symbolic_obstruction = bool(sp.simplify(tracefree_post_spatial_trace + 3 * q_shear) == 0 and q_positive_definite)
    all_samples_obstruct = all(row["energy_neutral_full_cancellation_impossible"] for row in rows)
    bounded_obstruction_passed = bool(p1974_present and p1977_no_go and energy_symbolic_obstruction and tracefree_symbolic_obstruction and all_samples_obstruct)

    out = {
        "ledger_id": "P1978_S928_STRICT_ENERGY_NEUTRAL_TENSOR_TRANSPORT_OBSTRUCTION",
        "packet_id": "P1978",
        "stage_id": "S928",
        "produced_by": "p1978_s928_strict_energy_neutral_tensor_transport_obstruction.py",
        "timestamp_utc": DEFAULT_TIMESTAMP_UTC,
        "route": "strict_only",
        "depends_on": {
            "p1974_present": p1974_present,
            "p1977_bounded_no_go_present": p1977_no_go,
        },
        "bounded_assumptions": {
            "A1_residual_source": "Use the P1974 diagonal Bianchi-I residual vector in components (G00,G11/a1^2,G22/a2^2,G33/a3^2).",
            "A2_energy_neutral_transport": "The non-minimal tensorial transport is spatial/shear-only and has no lapse/energy component: U00=0.",
            "A3_tracefree_subcase": "For the stricter shear-only subcase, the spatial transport is trace-free: U11+U22+U33=0.",
            "A4_nonzero_trace_free_shear": "(sigma1,sigma2) != (0,0), so Q_shear>0 by the positive-definite quadratic form.",
        },
        "symbolic_core": {
            "sigma3": str(sigma3),
            "dsigma3": str(dsigma3),
            "Q_shear": str(q_shear),
            "S_traceless_components": [str(s1), str(s2), str(s3)],
            "S_traceless_sum": str(sp.factor(s1 + s2 + s3)),
            "P1974_residual_vector": [str(expr) for expr in residual],
            "required_U00_for_full_cancellation": str(required_u0_for_full_cancellation),
            "energy_neutral_condition": str(energy_neutral_condition),
            "energy_neutral_leftover_G00": str(energy_neutral_energy_residual),
            "energy_neutral_full_cancellation_requires_Q_shear_zero": str(energy_neutral_requires_q_zero),
            "spatial_trace_residual": str(spatial_trace_residual),
            "post_spatial_trace_general": str(post_spatial_trace),
            "tracefree_spatial_leftover_sum": str(tracefree_post_spatial_trace),
            "Q_matrix": str(q_matrix),
            "Q_eigenvalues_exact": ["1/2", "3/2"],
            "Q_eigenvalues_scipy": [float(ev) for ev in scipy_eigs],
            "Q_positive_definite": q_positive_definite,
        },
        "numeric_replay_table": rows,
        "gatekeeper_checks": {
            "energy_neutral_leaves_G00_minus_Q": energy_symbolic_obstruction,
            "tracefree_spatial_leaves_sum_minus_3Q": tracefree_symbolic_obstruction,
            "q_shear_positive_definite": q_positive_definite,
            "all_nonzero_samples_obstruct_energy_neutral_cancellation": all_samples_obstruct,
            "bounded_obstruction_passed": bounded_obstruction_passed,
        },
        "result_kind": "PASS_BOUNDED_OBSTRUCTION" if bounded_obstruction_passed else "OPEN_OBSTRUCTION_WITH_TRACE",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "obstruction_tags": [
            "ENERGY_NEUTRAL_TENSOR_TRANSPORT_CANNOT_CANCEL_BIANCHI_I_G00_RESIDUAL",
            "TRACEFREE_SPATIAL_TRANSPORT_CANNOT_CANCEL_SPATIAL_TRACE_RESIDUAL",
            "BACKGROUND_INDEPENDENCE_REQUIRES_LAPSE_ENERGY_OR_DERIVED_GRAVITATIONAL_SHEAR_SECTOR",
            "GLOBAL_TOE_CLOSURE_STILL_OPEN",
        ],
        "theorem_export": {
            "bounded_negative_statement": "For the P1974 residual, any non-minimal spatial/shear-only transport with U00=0 leaves the G00 residual -Q_shear.  Since Q_shear is positive definite for nonzero trace-free shear, full componentwise FRW/Bianchi-I cancellation is impossible in this energy-neutral transport class.",
            "tracefree_subcase_statement": "If the spatial transport is also trace-free, the spatial trace residual remains -3*Q_shear, giving an independent obstruction.",
            "remaining_escape_routes": [
                "derive a strict lapse/energy component from a gravitational/shear variational sector and rerun the energy-sign audit",
                "derive a non-energy-neutral tensorial connection whose U00 contribution is theorem-grade and not an ad hoc negative-energy source",
                "prove that the current strict variational operator bundle cannot export such a lapse/energy contribution",
            ],
            "not_licensed": [
                "global no-go theorem for all tensorial transports",
                "negative-energy admission",
                "background-independence closure",
                "PO2/PO3 closure",
                "BRST or Cutkosky closure",
                "QW-2191 selector discharge",
                "ToE closure",
            ],
        },
        "false_pass_guard": "P1978 only rules out energy-neutral spatial/shear-only non-minimal transport for the P1974 residual.  It must not be read as a global no-go theorem or as background-independence closure.",
        "next_honest_step": "Test the remaining non-energy-neutral escape route: construct a strict variational shear/lapse provider with U00=-Q_shear, or prove that the current strict variational operator bundle cannot export that lapse/energy term.",
        "lay_explanation": "Sprawdziliśmy, czy da się naprawić błąd Bianchi-I samym geometrycznym przeniesieniem w częściach przestrzennych, bez dokładania energii.  Nie da się: równanie energii G00 zostawia składnik -Q_shear.  To zawęża drogę ToE do trudniejszego zadania: trzeba wyprowadzić prawdziwy sektor geometryczny z wkładem energii/lapse, a nie tylko przesunąć indeksy przestrzenne.",
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
