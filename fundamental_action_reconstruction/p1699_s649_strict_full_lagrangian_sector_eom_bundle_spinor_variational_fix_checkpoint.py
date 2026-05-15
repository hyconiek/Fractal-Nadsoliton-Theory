#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1698 = GEN / "p1698_s648_strict_full_lagrangian_sector_eom_bundle_nonproxy_scaffold_checkpoint.json"
OUT = GEN / "p1699_s649_strict_full_lagrangian_sector_eom_bundle_spinor_variational_fix_checkpoint.json"


def main() -> None:
    p1698 = json.loads(IN1698.read_text(encoding="utf-8"))
    coeff = p1698.get("forward_coefficient_map_symbolic", {})

    x = sp.symbols("x", real=True)
    A = sp.Function("A")(x)
    h = sp.Function("h")(x)
    phi = sp.Function("phi")(x)
    psi = sp.Function("psi")(x)
    psib = sp.Function("psib")(x)

    R0 = sp.symbols("R0", real=True)

    Z1 = sp.Float(coeff.get("Z1", "1.0"))
    muH2 = sp.Float(coeff.get("muH2", "1.0"))
    lambdaH = sp.Float(coeff.get("lambdaH", "1.0"))
    xiHR = sp.Float(coeff.get("xiHR", "0.0"))
    cR2 = sp.Float(coeff.get("cR2", "0.0"))
    mpsi = sp.Float("1.0")
    ypsi = sp.Float("1.0")

    L_gauge = sp.expand(Z1 * sp.diff(A, x) ** 2 / 2)
    L_higgs = sp.expand(
        sp.diff(h, x) ** 2 / 2
        - muH2 * h**2 / 2
        - lambdaH * h**4 / 4
        - xiHR * R0 * h**2 / 2
    )
    L_scalar = sp.expand(sp.diff(phi, x) ** 2 / 2 - cR2 * phi**4 / 4)

    # fermion proxy fixed: use independent psi and psib fields + hermitian-symmetrized kinetic term
    L_fermion = sp.expand(
        sp.I * (psib * sp.diff(psi, x) - sp.diff(psib, x) * psi) / 2
        - (mpsi + ypsi * h) * psib * psi
    )
    L_mix = sp.expand(-h**2 * phi**2 / 2)

    L_total_reduced = sp.expand(L_gauge + L_higgs + L_scalar + L_fermion + L_mix)

    def el(L: sp.Expr, q: sp.Expr) -> sp.Expr:
        return sp.simplify(sp.diff(sp.diff(L, sp.diff(q, x)), x) - sp.diff(L, q))

    eom_A = el(L_total_reduced, A)
    eom_h = el(L_total_reduced, h)
    eom_phi = el(L_total_reduced, phi)
    eom_psi = el(L_total_reduced, psi)
    eom_psib = el(L_total_reduced, psib)

    payload = {
        "checkpoint": "P1699_S649",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "source_checkpoint": "P1698_S648",
        "chain": "K_strict -> coefficients -> full explicit L_total anchor -> corrected spinor variational reduced bundle -> EOM",
        "forward_coefficient_map_symbolic": coeff,
        "full_lagrangian_density_explicit_anchor": p1698.get("full_lagrangian_density_explicit_anchor", {}),
        "reduced_bundle_with_spinor_fix": {
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
                "EOM_psib": str(eom_psib),
            },
        },
        "consistency_note": "Fermion EL no longer collapses to zero artifact from single-field total-derivative proxy; psi/psib pair tracks Dirac-like variational structure in reduced form.",
        "reverse_direction_status": {
            "eom_to_variational_origin": "OPEN_THEOREM_REQUIRED_NONPROXY",
            "helmholtz_link": "PARTIAL_REDUCED_PROXY_ONLY",
        },
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "open_obligations": [
            "full_covariant_spinor_bundle_nonproxy",
            "full_metric_tensor_nonproxy",
            "counterterm_flow_closure_theorem",
            "BRST_nilpotency_and_cohomology_theorem",
            "Cutkosky_unitarity_theorem_full_sector",
            "background_family_independence_theorem",
            "strict_selector_source_or_symmetry_breaking_premise_for_QW2191",
        ],
        "next_honest_step": "Zastąpić reduced spinor proxy pełnym kowariantnym spinor bundle z połączeniem spinowym i gamma-macierzami, a następnie domknąć theorem-level QG closures.",
        "lay_summary": "Naprawiliśmy techniczny problem w sektorze fermionowym: wcześniejsza uproszczona forma dawała sztucznie zerowe równanie. Teraz równania fermionowe są niezerowe i bardziej fizycznie wiarygodne, ale to nadal etap pośredni.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
