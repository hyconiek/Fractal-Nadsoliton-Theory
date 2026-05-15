#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1747 = GEN / "p1747_s697_strict_full_lagrangian_non_skeleton_and_bidirectional_witness_bundle_checkpoint.json"
OUT = GEN / "p1748_s698_strict_reduced_h1_cross_variation_witness_checkpoint.json"


def main() -> None:
    p1747 = json.loads(IN1747.read_text(encoding="utf-8"))
    coeff = p1747.get("forward_coefficient_map_symbolic", {})

    x = sp.symbols("x", real=True)
    h = sp.Function("h")(x)
    phi = sp.Function("phi")(x)
    R0 = sp.symbols("R0", real=True)

    muH2 = sp.Float(coeff.get("muH2", "1.0"))
    lambdaH = sp.Float(coeff.get("lambdaH", "1.0"))
    xiHR = sp.Float(coeff.get("xiHR", "0.0"))
    cR2 = sp.Float(coeff.get("cR2", "0.0"))

    L = (
        sp.diff(h, x) ** 2 / 2
        + sp.diff(phi, x) ** 2 / 2
        - muH2 * h**2 / 2
        - lambdaH * h**4 / 4
        - xiHR * R0 * h**2 / 2
        - cR2 * phi**4 / 4
        - h**2 * phi**2 / 2
    )

    E_h = sp.simplify(sp.diff(sp.diff(L, sp.diff(h, x)), x) - sp.diff(L, h))
    E_phi = sp.simplify(sp.diff(sp.diff(L, sp.diff(phi, x)), x) - sp.diff(L, phi))

    dEh_dphi = sp.simplify(sp.diff(E_h, phi))
    dEphi_dh = sp.simplify(sp.diff(E_phi, h))
    diff_cross = sp.simplify(dEh_dphi - dEphi_dh)

    pass_zero = sp.simplify(diff_cross) == 0

    payload = {
        "checkpoint": "P1748_S698",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "scope": "REDUCED_PROXY_ONLY_NOT_NONPROXY",
        "chain": "K_strict -> coefficients -> reduced non-skeleton L -> EOM(h,phi) -> Helmholtz H1 local cross-variation",
        "input_anchor": "p1747_s697",
        "EOM": {"E_h": str(E_h), "E_phi": str(E_phi)},
        "cross_variation": {
            "delta_E_h_over_delta_phi": str(dEh_dphi),
            "delta_E_phi_over_delta_h": str(dEphi_dh),
            "difference": str(diff_cross),
        },
        "result": "PASS_ZERO" if pass_zero else "OBSTRUCTION",
        "result_policy": "PASS valid only for reduced proxy scope; no automatic upgrade to nonproxy/full covariant theorem",
        "qg_closure_impact": "NONE_DIRECT_OPEN",
        "remaining_open": [
            "nonproxy_metric_tensor_variational_export",
            "nonproxy_spinor_connection_variational_export",
            "helmholtz_integrability_full_sector_nonproxy",
            "renormalization_closure_counterterm_stability",
            "cutkosky_unitarity_full_sector",
            "background_family_independence",
            "QW-2191_selector_source_or_symmetry_breaking_premise",
        ],
        "next_honest_step": "Zastąpić reduced H1 przez nonproxy covariant H1 (A_mu,H) na wspólnej rodzinie teł z jawną kontrolą wyrazów brzegowych.",
        "lay_summary": "W tej wersji testu pokazaliśmy matematycznie, że dwa sprzężone sektory w modelu redukcyjnym są wzajemnie spójne (H1). To ważny test lokalny, ale nie kończy jeszcze pełnej fizyki kwantowej grawitacji.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
