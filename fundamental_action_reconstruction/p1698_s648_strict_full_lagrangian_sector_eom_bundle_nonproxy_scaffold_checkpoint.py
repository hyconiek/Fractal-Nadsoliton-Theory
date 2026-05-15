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
OUT = GEN / "p1698_s648_strict_full_lagrangian_sector_eom_bundle_nonproxy_scaffold_checkpoint.json"


def _f(coeff: dict[str, str], key: str, default: str) -> sp.Float:
    return sp.Float(coeff.get(key, default))


def main() -> None:
    p1694 = json.loads(IN1694.read_text(encoding="utf-8"))
    p1662 = json.loads(IN1662.read_text(encoding="utf-8"))
    coeff = p1694.get("forward_coefficient_map_symbolic", {})

    x = sp.symbols("x", real=True)
    A = sp.Function("A")(x)
    h = sp.Function("h")(x)
    phi = sp.Function("phi")(x)
    psi = sp.Function("psi")(x)

    R0 = sp.symbols("R0", real=True)

    Z1 = _f(coeff, "Z1", "1.0")
    muH2 = _f(coeff, "muH2", "1.0")
    lambdaH = _f(coeff, "lambdaH", "1.0")
    xiHR = _f(coeff, "xiHR", "0.0")
    cR2 = _f(coeff, "cR2", "0.0")

    # strict-only reduced 1D proxy for full-sector EL export scaffold
    # (explicitly not a replacement for full tensor/spinor nonproxy objects)
    L_gauge = sp.expand(Z1 * sp.diff(A, x) ** 2 / 2)
    L_higgs = sp.expand(
        sp.diff(h, x) ** 2 / 2
        - muH2 * h**2 / 2
        - lambdaH * h**4 / 4
        - xiHR * R0 * h**2 / 2
    )
    L_scalar = sp.expand(sp.diff(phi, x) ** 2 / 2 - cR2 * phi**4 / 4)
    L_fermion = sp.expand(sp.I * psi * sp.diff(psi, x))
    L_mix = sp.expand(-h**2 * phi**2 / 2)

    L_total_reduced = sp.expand(L_gauge + L_higgs + L_scalar + L_fermion + L_mix)

    eom_A = sp.simplify(sp.diff(sp.diff(L_total_reduced, sp.diff(A, x)), x) - sp.diff(L_total_reduced, A))
    eom_h = sp.simplify(sp.diff(sp.diff(L_total_reduced, sp.diff(h, x)), x) - sp.diff(L_total_reduced, h))
    eom_phi = sp.simplify(sp.diff(sp.diff(L_total_reduced, sp.diff(phi, x)), x) - sp.diff(L_total_reduced, phi))
    eom_psi = sp.simplify(sp.diff(sp.diff(L_total_reduced, sp.diff(psi, x)), x) - sp.diff(L_total_reduced, psi))

    payload = {
        "checkpoint": "P1698_S648",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "chain": "K_strict -> coefficients -> full explicit L_total anchor -> sector EOM bundle scaffold",
        "kernel": p1694.get("kernel", "K_strict(d)=cos(omega*d+phi)/(1+beta*d^eta)"),
        "kernel_input": p1694.get("kernel_input", {}),
        "forward_coefficient_map_symbolic": coeff,
        "full_lagrangian_density_explicit_anchor": p1662.get("full_lagrangian_density_explicit", {}),
        "reduced_bundle_for_symbolic_el_scaffold": {
            "background_marker": "R0",
            "L_total_reduced": str(L_total_reduced),
            "L_sector": {
                "L_gauge": str(L_gauge),
                "L_higgs": str(L_higgs),
                "L_scalar": str(L_scalar),
                "L_fermion": str(L_fermion),
                "L_mix": str(L_mix),
            },
            "EOM": {
                "EOM_A": str(eom_A),
                "EOM_h": str(eom_h),
                "EOM_phi": str(eom_phi),
                "EOM_psi": str(eom_psi),
            },
        },
        "reverse_direction_status": {
            "eom_to_variational_origin": "OPEN_THEOREM_REQUIRED",
            "helmholtz_witness_link": "REQUIRES_NONPROXY_TENSOR_SPINOR_EXPORT",
        },
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "open_obligations": [
            "full_covariant_tensor_export_metric_nonproxy",
            "full_covariant_spinor_export_gamma_connection_nonproxy",
            "counterterm_flow_closure_theorem",
            "BRST_nilpotency_and_cohomology_theorem",
            "Cutkosky_unitarity_theorem_full_sector",
            "background_family_independence_theorem",
            "strict_selector_source_or_symmetry_breaking_premise_for_QW2191",
        ],
        "next_honest_step": "Przejść z reduced scaffold do pełnego nieproxy bundle tensor+spinor, a następnie zamknąć theorem-level tor QG (renormalizacja, unitarność, background independence) na jednej rodzinie teł.",
        "lay_summary": "Wyprowadziliśmy spójny pakiet równań ruchu z pełnego strict łańcucha, ale na poziomie scaffold. To przybliża pełny model, lecz finalny certyfikat wymaga jeszcze pełnych obiektów kowariantnych i dowodów kwantowej grawitacji.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
