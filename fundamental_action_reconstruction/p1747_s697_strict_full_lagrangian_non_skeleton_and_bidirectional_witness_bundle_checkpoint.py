#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1694 = GEN / "p1694_s644_strict_kernel_to_full_lagrangian_bidirectional_map_witness.json"
IN1662 = GEN / "p1662_s612_strict_full_lagrangian_explicit_density_summary.json"
OUT = GEN / "p1747_s697_strict_full_lagrangian_non_skeleton_and_bidirectional_witness_bundle_checkpoint.json"


def _f(coeff: dict[str, str], key: str, default: str) -> sp.Float:
    return sp.Float(coeff.get(key, default))


def main() -> None:
    p1694 = json.loads(IN1694.read_text(encoding="utf-8"))
    p1662 = json.loads(IN1662.read_text(encoding="utf-8"))
    coeff = p1694.get("forward_coefficient_map_symbolic", {})

    x = sp.symbols("x", real=True)
    A, h, phi = [sp.Function(n)(x) for n in ("A", "h", "phi")]
    psi = sp.Function("psi")(x)
    psib = sp.Function("psib")(x)

    R0 = sp.symbols("R0", real=True)

    Z1 = _f(coeff, "Z1", "1.0")
    muH2 = _f(coeff, "muH2", "1.0")
    lambdaH = _f(coeff, "lambdaH", "1.0")
    xiHR = _f(coeff, "xiHR", "0.0")
    cR2 = _f(coeff, "cR2", "0.0")

    # Non-skeleton strict reduced proxy with independent psi/psib (no false zero artifact)
    L_gauge = sp.expand(Z1 * sp.diff(A, x) ** 2 / 2)
    L_higgs = sp.expand(sp.diff(h, x) ** 2 / 2 - muH2 * h**2 / 2 - lambdaH * h**4 / 4 - xiHR * R0 * h**2 / 2)
    L_scalar = sp.expand(sp.diff(phi, x) ** 2 / 2 - cR2 * phi**4 / 4)
    L_mix = sp.expand(-h**2 * phi**2 / 2)
    L_fermion = sp.expand(sp.I / 2 * (psib * sp.diff(psi, x) - sp.diff(psib, x) * psi))

    L_total = sp.expand(L_gauge + L_higgs + L_scalar + L_mix + L_fermion)

    def eom(field: sp.Expr) -> sp.Expr:
        return sp.simplify(sp.diff(sp.diff(L_total, sp.diff(field, x)), x) - sp.diff(L_total, field))

    eoms = {
        "EOM_A": str(eom(A)),
        "EOM_h": str(eom(h)),
        "EOM_phi": str(eom(phi)),
        "EOM_psi": str(eom(psi)),
        "EOM_psib": str(eom(psib)),
    }

    payload = {
        "checkpoint": "P1747_S697",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "chain": "K_strict -> coefficients -> non-skeleton L_total export -> EOM bundle -> reverse witness obligations",
        "kernel": p1694.get("kernel", "K_strict(d)=cos(omega*d+phi)/(1+beta*d^eta)"),
        "kernel_input": p1694.get("kernel_input", {}),
        "forward_coefficient_map_symbolic": coeff,
        "full_lagrangian_density_explicit_anchor": p1662.get("full_lagrangian_density_explicit", {}),
        "strict_non_skeleton_reduced_proxy": {
            "background_marker": "R0",
            "L_total_reduced": str(L_total),
            "L_sector": {
                "L_gauge": str(L_gauge),
                "L_higgs": str(L_higgs),
                "L_scalar": str(L_scalar),
                "L_mix": str(L_mix),
                "L_fermion_hermitian_independent_psi_psib": str(L_fermion),
            },
            "EOM_bundle": eoms,
        },
        "strict_core_closure_gate": {
            "direction_forward_kernel_to_eom": "EXPORTED_REDUCED_NON_SKELETON_PROXY",
            "direction_reverse_eom_to_action": "OPEN_WITNESS_CHAIN_REQUIRED",
            "required_missing_exports": [
                "nonproxy_metric_tensor_variational_export",
                "nonproxy_spinor_connection_variational_export",
                "nonproxy_BRST_cohomology_export",
            ],
            "required_missing_theorems": [
                "helmholtz_integrability_full_sector",
                "renormalization_closure_counterterm_stability",
                "cutkosky_unitarity_full_sector",
                "background_family_independence",
            ],
            "selector_boundary": "QW-2191 remains OPEN unless explicit strict selector source/symmetry-breaking premise is exported",
            "final_status": "NO_FALSE_PASS_STRICT_CORE_CLOSURE_NOT_YET_DISCHARGED",
        },
        "next_honest_step": "Dostarczyć pierwszy nieproxy witness reverse: Helmholtz H1 dla sektora gauge+scalar z jawnym eksportem operatora Jacobian symmetry test.",
        "lay_summary": "Mamy pełniejszy (nieszkieletowy) ścisły Lagrangian redukcyjny i równania ruchu, ale pełne domknięcie ToE nadal wymaga brakujących eksportów i twierdzeń reverse (renormalizacja, unitarność, niezależność od tła).",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
