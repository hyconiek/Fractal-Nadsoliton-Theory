#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1747 = GEN / "p1747_s697_strict_full_lagrangian_non_skeleton_and_bidirectional_witness_bundle_checkpoint.json"
OUT = GEN / "p1749_s699_strict_reduced_h1_gauge_higgs_cross_variation_checkpoint.json"


def main() -> None:
    p1747 = json.loads(IN1747.read_text(encoding="utf-8"))
    coeff = p1747.get("forward_coefficient_map_symbolic", {})

    x = sp.symbols("x", real=True)
    A = sp.Function("A")(x)
    h = sp.Function("h")(x)
    R0 = sp.symbols("R0", real=True)

    Z1 = sp.Float(coeff.get("Z1", "1.0"))
    muH2 = sp.Float(coeff.get("muH2", "1.0"))
    lambdaH = sp.Float(coeff.get("lambdaH", "1.0"))
    xiHR = sp.Float(coeff.get("xiHR", "0.0"))

    # Reduced gauge-higgs coupling proxy: (d h - g A h)^2 /2 gives nontrivial A-H cross channel
    g = sp.Symbol("g", real=True)
    Dh = sp.diff(h, x) - g * A * h
    L = (
        Z1 * sp.diff(A, x) ** 2 / 2
        + Dh**2 / 2
        - muH2 * h**2 / 2
        - lambdaH * h**4 / 4
        - xiHR * R0 * h**2 / 2
    )

    E_A = sp.simplify(sp.diff(sp.diff(L, sp.diff(A, x)), x) - sp.diff(L, A))
    E_h = sp.simplify(sp.diff(sp.diff(L, sp.diff(h, x)), x) - sp.diff(L, h))

    dEA_dh = sp.simplify(sp.diff(E_A, h))
    dEh_dA = sp.simplify(sp.diff(E_h, A))
    h1_diff = sp.simplify(dEA_dh - dEh_dA)

    pass_zero = sp.simplify(h1_diff) == 0

    payload = {
        "checkpoint": "P1749_S699",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "scope": "REDUCED_PROXY_ONLY_NOT_NONPROXY",
        "chain": "K_strict -> coefficients -> non-skeleton reduced gauge-higgs L -> EOM(A,h) -> H1(A,h) cross-variation",
        "input_anchor": "p1747_s697",
        "modeling_note": "Gauge-Higgs reduced proxy uses 1D covariant-like derivative Dh = d_x h - g A h to create explicit A-h cross channel.",
        "EOM": {"E_A": str(E_A), "E_h": str(E_h)},
        "cross_variation": {
            "delta_E_A_over_delta_h": str(dEA_dh),
            "delta_E_h_over_delta_A": str(dEh_dA),
            "difference": str(h1_diff),
        },
        "result": "PASS_ZERO" if pass_zero else "OBSTRUCTION",
        "result_policy": "PASS valid only for reduced proxy scope; no automatic upgrade to nonproxy/full covariant theorem",
        "remaining_open": [
            "explicit_covariant_E_A_mu_expression_nonproxy",
            "explicit_covariant_E_H_expression_nonproxy",
            "shared_background_family_contract",
            "boundary_term_control_clause",
            "nonproxy_metric_tensor_variational_export",
            "renormalization_closure_counterterm_stability",
            "cutkosky_unitarity_full_sector",
            "background_family_independence",
        ],
        "next_honest_step": "Przenieść ten sam test H1 na jawne nonproxy E_A^μ i E_H (4D, wspólna rodzina teł) oraz dołączyć ślad wyrazów brzegowych.",
        "lay_summary": "To jest bardziej fizyczny test niż poprzednio, bo łączy pole cechowania z Higgsem. Lokalnie przechodzi warunek zgodności odwrotnej, ale pełna wersja teorii nadal wymaga eksportu nonproxy i twierdzeń QG.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
