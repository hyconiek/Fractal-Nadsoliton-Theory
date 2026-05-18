#!/usr/bin/env python3
"""P1980 S930 strict ADM/Bianchi-I Einstein-Hilbert lapse variation witness.

This is the first non-token ADM calculation after P1979.  It keeps the lapse
N(t) in a diagonal Bianchi-I minisuperspace reduction of the Einstein-Hilbert
ADM bulk term and varies with respect to N.  The calculation shows that the
EH lapse equation already contains the geometric shear-energy contribution
G_nn(Bianchi-I)-G_nn(FRW) = -Q_shear.

Scope: EH/lapse equation only.  It is not the full curvature-squared variation,
not spatial EOM transport, and not background-independence closure.
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
OUT = GEN / "p1980_s930_strict_adm_bianchi_eh_lapse_variation_witness.json"


def load(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p1979 = load("p1979_s929_strict_lapse_shear_provider_export_audit.json")

    N, V = sp.symbols("N V", positive=True, real=True)
    H, sigma1, sigma2 = sp.symbols("H sigma1 sigma2", real=True)
    v1, v2, v3 = sp.symbols("v1 v2 v3", real=True)

    # ADM diagonal Bianchi-I with zero spatial curvature, no boundary term:
    # L_EH_bulk / prefactor = N*V*(K_ij K^ij - K^2).
    # Here v_i = dot(a_i)/a_i are coordinate-time logarithmic rates and
    # K^i_i = v_i/N.
    k_diag = [v1 / N, v2 / N, v3 / N]
    k_trace = sp.factor(sum(k_diag))
    kij_kij = sp.factor(sum(k**2 for k in k_diag))
    adm_bulk_lagrangian = sp.factor(N * V * (kij_kij - k_trace**2))
    lapse_variation = sp.factor(sp.diff(adm_bulk_lagrangian, N))
    # Remove the positive conventional factor 2V/N^2 to get the normalized
    # Hamiltonian/lapse constraint sum_{i<j} v_i v_j.
    normalized_lapse_constraint_v = sp.factor(lapse_variation * N**2 / (2 * V))
    expected_pair_sum_v = sp.factor(v1 * v2 + v1 * v3 + v2 * v3)
    variation_matches_pair_sum = sp.simplify(normalized_lapse_constraint_v - expected_pair_sum_v) == 0

    sigma3 = sp.factor(-sigma1 - sigma2)
    h1 = H + sigma1
    h2 = H + sigma2
    h3 = H + sigma3
    pair_sum_bianchi = sp.factor((h1 * h2 + h1 * h3 + h2 * h3))
    pair_sum_frw = sp.factor(3 * H**2)
    q_shear = sp.factor(sigma1**2 + sigma1 * sigma2 + sigma2**2)
    bianchi_minus_frw = sp.factor(pair_sum_bianchi - pair_sum_frw)
    residual_matches_p1974_g00 = sp.simplify(bianchi_minus_frw + q_shear) == 0

    q_matrix = sp.Matrix([[1, sp.Rational(1, 2)], [sp.Rational(1, 2), 1]])
    scipy_eigs = la.eigvalsh(np.array(q_matrix.tolist(), dtype=float))
    q_positive_definite = all(float(ev) > 0 for ev in scipy_eigs)

    sample_points = [
        {H: sp.Rational(1, 1), sigma1: sp.Rational(1, 10), sigma2: sp.Rational(-1, 20)},
        {H: sp.Rational(3, 2), sigma1: sp.Rational(1, 5), sigma2: sp.Rational(1, 10)},
        {H: sp.Rational(4, 5), sigma1: sp.Rational(-1, 8), sigma2: sp.Rational(1, 16)},
        {H: sp.Rational(5, 4), sigma1: sp.Rational(2, 7), sigma2: sp.Rational(-1, 9)},
    ]
    rows: list[dict[str, Any]] = []
    for idx, point in enumerate(sample_points):
        q_val = sp.simplify(q_shear.subs(point))
        bi_val = sp.simplify(pair_sum_bianchi.subs(point))
        frw_val = sp.simplify(pair_sum_frw.subs(point))
        res_val = sp.simplify(bianchi_minus_frw.subs(point))
        rows.append(
            {
                "sample_id": f"eh_lapse_variation_sample_{idx}",
                "input": {str(k): str(v) for k, v in point.items()},
                "Q_shear": str(q_val),
                "Gnn_BianchiI_EH_constraint": str(bi_val),
                "Gnn_FRW_EH_constraint": str(frw_val),
                "Bianchi_minus_FRW": str(res_val),
                "equals_minus_Q_shear": bool(sp.simplify(res_val + q_val) == 0),
            }
        )

    p1979_pass = p1979.get("result_kind") == "PASS_CURRENT_EXPORT_NONAVAILABILITY_AUDIT"
    all_samples_match = all(row["equals_minus_Q_shear"] for row in rows)
    theorem_pass = bool(p1979_pass and variation_matches_pair_sum and residual_matches_p1974_g00 and q_positive_definite and all_samples_match)

    out = {
        "ledger_id": "P1980_S930_STRICT_ADM_BIANCHI_EH_LAPSE_VARIATION_WITNESS",
        "packet_id": "P1980",
        "stage_id": "S930",
        "produced_by": "p1980_s930_strict_adm_bianchi_eh_lapse_variation_witness.py",
        "timestamp_utc": DEFAULT_TIMESTAMP_UTC,
        "route": "strict_only",
        "depends_on": {
            "p1979_current_export_audit_present": p1979_pass,
        },
        "adm_setup": {
            "metric": "ds^2 = -N(t)^2 dt^2 + a1(t)^2 dx^2 + a2(t)^2 dy^2 + a3(t)^2 dz^2",
            "spatial_curvature_R3": "0 for diagonal Bianchi-I",
            "K_i^i": ["v1/N", "v2/N", "v3/N"],
            "EH_bulk_density_prefactor_removed": "N*V*(K_ij*K^ij - K^2)",
            "boundary_terms_status": "not needed for lapse algebra of the bulk ADM constraint in this minisuperspace witness",
        },
        "symbolic_variation": {
            "K_trace": str(k_trace),
            "KijKij": str(kij_kij),
            "ADM_EH_bulk_lagrangian": str(adm_bulk_lagrangian),
            "dL_dN": str(lapse_variation),
            "normalized_lapse_constraint_v": str(normalized_lapse_constraint_v),
            "expected_pair_sum_v": str(expected_pair_sum_v),
            "variation_matches_pair_sum": variation_matches_pair_sum,
        },
        "shear_specialization": {
            "sigma3": str(sigma3),
            "H_i": [str(h1), str(h2), str(h3)],
            "Q_shear": str(q_shear),
            "Gnn_BianchiI_EH_constraint": str(pair_sum_bianchi),
            "Gnn_FRW_EH_constraint": str(pair_sum_frw),
            "Bianchi_minus_FRW": str(bianchi_minus_frw),
            "matches_P1974_G00_residual_minus_Q_shear": residual_matches_p1974_g00,
            "Q_matrix": str(q_matrix),
            "Q_eigenvalues_exact": ["1/2", "3/2"],
            "Q_eigenvalues_scipy": [float(ev) for ev in scipy_eigs],
            "Q_positive_definite": q_positive_definite,
        },
        "numeric_replay_table": rows,
        "gatekeeper_checks": {
            "p1979_identified_missing_export_before_variation": p1979_pass,
            "adm_lapse_variation_matches_pair_sum": variation_matches_pair_sum,
            "eh_bianchi_minus_frw_equals_minus_q_shear": residual_matches_p1974_g00,
            "q_shear_positive_definite": q_positive_definite,
            "all_numeric_replay_samples_match_minus_q": all_samples_match,
            "eh_lapse_witness_passed": theorem_pass,
        },
        "result_kind": "PASS_EH_LAPSE_SHEAR_TERM_DERIVED" if theorem_pass else "OPEN_OBSTRUCTION_WITH_TRACE",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "obstruction_tags": [
            "EH_LAPSE_EQUATION_DERIVES_GEOMETRIC_MINUS_Q_SHEAR_TERM",
            "CURVATURE_SQUARED_ADM_VARIATION_STILL_MISSING",
            "SPATIAL_EOM_TRANSPORT_STILL_MISSING",
            "BACKGROUND_INDEPENDENCE_REMAINS_OPEN",
            "GLOBAL_TOE_CLOSURE_STILL_OPEN",
        ],
        "theorem_export": {
            "positive_statement": "For diagonal Bianchi-I in the ADM Einstein-Hilbert bulk term, lapse variation gives the normalized constraint H1*H2+H1*H3+H2*H3.  Under trace-free shear H_i=H+sigma_i, this equals 3H^2-Q_shear, so the Bianchi-I minus FRW lapse residual is exactly -Q_shear.",
            "correction_to_P1979_audit_reading": "P1979 was only a current-export token/certificate audit.  P1980 shows that the EH covariant term does carry the required geometric shear-energy sign after explicit ADM lapse variation, even though the registry did not export that certificate textually.",
            "not_licensed": [
                "full curvature-squared ADM variation",
                "spatial Bianchi-I EOM cancellation",
                "global background-independence closure",
                "PO2/PO3 closure",
                "BRST or Cutkosky closure",
                "QW-2191 selector discharge",
                "ToE closure",
            ],
        },
        "false_pass_guard": "P1980 is a positive EH lapse-equation witness only.  It must not be used to claim full background-independence, curvature-squared closure, spatial EOM closure, or ToE closure.",
        "next_honest_step": "Extend the ADM/Bianchi-I variation from the Einstein-Hilbert bulk term to the strict curvature-squared terms R^2, Ricci^2, Riemann^2, and G_GB, with lapse N(t) retained until variation, then compare their shear contributions against the remaining P1974/P1975 component obligations.",
        "lay_explanation": "Zamiast tylko szukać słów w rejestrze, wykonaliśmy pierwszy prawdziwy rachunek z funkcją lapse N(t).  Dla zwykłego składnika Einsteina-Hilberta anizotropia sama daje dokładnie brakujący znak -Q_shear w równaniu energii.  To dobry znak, ale jeszcze nie pełny dowód: trzeba policzyć trudniejsze składniki kwadratowe w krzywiźnie i równania przestrzenne.",
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
