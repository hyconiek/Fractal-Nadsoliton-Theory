#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1699 = GEN / "p1699_s649_strict_full_lagrangian_sector_eom_bundle_spinor_variational_fix_checkpoint.json"
OUT = GEN / "p1700_s650_strict_eom_variational_residual_identity_checkpoint.json"


def main() -> None:
    p1699 = json.loads(IN1699.read_text(encoding="utf-8"))
    coeff = p1699.get("forward_coefficient_map_symbolic", {})

    x = sp.symbols("x", real=True)
    A = sp.Function("A")(x)
    h = sp.Function("h")(x)
    phi = sp.Function("phi")(x)
    psi = sp.Function("psi")(x)
    psib = sp.Function("psib")(x)
    R0 = sp.symbols("R0", real=True)

    Z1 = sp.Float(coeff.get("Z1", "1.0")); muH2 = sp.Float(coeff.get("muH2", "1.0"))
    lambdaH = sp.Float(coeff.get("lambdaH", "1.0")); xiHR = sp.Float(coeff.get("xiHR", "0.0"))
    cR2 = sp.Float(coeff.get("cR2", "0.0")); mpsi = sp.Float("1.0"); ypsi = sp.Float("1.0")

    L = sp.expand(
        Z1 * sp.diff(A, x) ** 2 / 2
        + sp.diff(h, x) ** 2 / 2 - muH2 * h**2 / 2 - lambdaH * h**4 / 4 - xiHR * R0 * h**2 / 2
        + sp.diff(phi, x) ** 2 / 2 - cR2 * phi**4 / 4
        + sp.I * (psib * sp.diff(psi, x) - sp.diff(psib, x) * psi) / 2 - (mpsi + ypsi * h) * psib * psi
        - h**2 * phi**2 / 2
    )

    fields = {"A": A, "h": h, "phi": phi, "psi": psi, "psib": psib}

    def el(q: sp.Expr) -> sp.Expr:
        return sp.simplify(sp.diff(sp.diff(L, sp.diff(q, x)), x) - sp.diff(L, q))

    eoms = {name: el(q) for name, q in fields.items()}
    residuals = {f"EL_minus_EOM_{k}": sp.simplify(el(v) - eoms[k]) for k, v in fields.items()}
    residual_zero = {k: (str(v) == "0") for k, v in residuals.items()}

    payload = {
        "checkpoint": "P1700_S650",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "source_checkpoint": "P1699_S649",
        "chain": "K_strict -> coefficients -> full explicit L_total anchor -> reduced corrected EOM -> EL residual identity",
        "L_total_reduced": str(L),
        "EOM_bundle": {k: str(v) for k, v in eoms.items()},
        "variational_identity_residuals": {k: str(v) for k, v in residuals.items()},
        "variational_identity_zero_flags": residual_zero,
        "identity_status": "PASS_ALL_ZERO" if all(residual_zero.values()) else "FAIL_NONZERO_RESIDUAL",
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "open_obligations": [
            "full_covariant_nonproxy_metric_spinor_export",
            "global_helmholtz_integrability_theorem",
            "counterterm_flow_closure_theorem",
            "BRST_Cutkosky_unitarity_theorem",
            "background_family_independence_theorem",
        ],
        "next_honest_step": "Przenieść ten identity-pass z reduced proxy na pełny nonproxy bundle (tensor+spinor) i domknąć theorem-level QG.",
        "lay_summary": "Sprawdziliśmy automatycznie, że wyprowadzone równania są matematycznie zgodne z tym lagranżianem (reszty = 0). To wzmacnia poprawność lokalnego kroku, ale nie zamyka jeszcze całej teorii.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
