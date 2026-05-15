#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1664 = GEN / "p1664_s614_strict_full_lagrangian_manifest_and_inversion.json"
IN1691 = GEN / "p1691_s641_strict_full_chain_lagrangian_to_qg_theorem_obligation_matrix.json"
OUT = GEN / "p1692_s642_strict_sympy_coefficient_to_eom_witness.json"


def main() -> None:
    p1664 = json.loads(IN1664.read_text(encoding="utf-8"))
    p1691 = json.loads(IN1691.read_text(encoding="utf-8"))
    coeff = p1664["coefficient_map"]

    # strict kernel parameters and derived strict coefficients
    omega = sp.Float(p1664["kernel_input"]["omega"])
    beta = sp.Float(p1664["kernel_input"]["beta"])
    eta = sp.Float(p1664["kernel_input"]["eta"])
    A = sp.Float(p1664["kernel_input"]["A"])

    x = sp.symbols("x", real=True)
    phi = sp.Function("phi")(x)
    H = sp.Function("H")(x)

    mphi2 = A * omega**2 * (1 + beta)
    lam3 = A * beta / (1 + eta)
    lam4 = A * beta * eta / (1 + eta)
    lam_phiH = sp.Float(coeff["xiHR"]) * sp.Float(coeff["lambdaH"])

    # reduced 1D sector density for symbolic EL replay witness
    L_scalar = sp.Rational(1, 2) * sp.diff(phi, x) ** 2 - sp.Rational(1, 2) * mphi2 * phi**2 - lam3 * phi**3 / sp.Integer(6) - lam4 * phi**4 / sp.Integer(24)
    L_hmix = sp.Rational(1, 2) * sp.diff(H, x) ** 2 - sp.Float(coeff["muH2"]) * H**2 / 2 - sp.Float(coeff["lambdaH"]) * H**4 / 4 - lam_phiH * phi**2 * H**2 / 2
    L_red = sp.expand(L_scalar + L_hmix)

    eom_phi = sp.simplify(sp.diff(sp.diff(L_red, sp.diff(phi, x)), x) - sp.diff(L_red, phi))
    eom_H = sp.simplify(sp.diff(sp.diff(L_red, sp.diff(H, x)), x) - sp.diff(L_red, H))

    payload = {
        "checkpoint": "P1692_S642",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "chain": "K_strict -> coefficients -> full L_total (anchor) -> reduced symbolic sector -> Euler-Lagrange witness",
        "kernel_input": p1664["kernel_input"],
        "coefficient_anchor": coeff,
        "full_lagrangian_density_anchor": p1691.get("full_lagrangian_density", {}),
        "reduced_density": str(L_red),
        "symbolic_eom": {
            "phi": str(eom_phi),
            "H": str(eom_H),
        },
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "open_obligations": [
            "full_covariant_EL_export_for_all_SM_GR_sectors",
            "theorem_level_counterterm_flow_closure",
            "theorem_level_BRST_cohomology_closure",
            "theorem_level_Cutkosky_unitarity_full_sector",
            "background_family_independence_theorem",
        ],
        "next_honest_step": "Podnieść replay SymPy z sektora zredukowanego do pełnego kowariantnego eksportu EL dla całego L_total i spiąć z theorem-level QG witnessami.",
        "lay_summary": "Policzyliśmy symbolicznie równania ruchu z jawnej części lagranżjanu, żeby łańcuch od kernela strict do EOM był mechanicznie sprawdzalny. To nadal nie zamyka problemu kwantowej grawitacji.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
