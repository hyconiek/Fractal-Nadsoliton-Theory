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
OUT = GEN / "p1694_s644_strict_kernel_to_full_lagrangian_bidirectional_map_witness.json"


def main() -> None:
    p1664 = json.loads(IN1664.read_text(encoding="utf-8"))
    p1691 = json.loads(IN1691.read_text(encoding="utf-8"))

    kin = p1664["kernel_input"]
    coeff = p1664["coefficient_map"]

    omega = sp.Float(kin["omega"])
    beta = sp.Float(kin["beta"])
    eta = sp.Float(kin["eta"])
    A = sp.Float(kin["A"])

    # forward: strict kernel -> coefficients (anchor formulas)
    fwd = {
        "Mpl2": sp.simplify(A * (1 + beta)),
        "cR2": sp.simplify(A * beta / (1 + eta)),
        "cRic2": sp.simplify(A * beta * eta / (1 + eta)),
        "cRiem2": sp.simplify(A * beta * (1 + eta) / 4),
        "muH2": sp.simplify(A * omega**2),
        "lambdaH": sp.simplify((1 + eta**2) / (1 + beta)),
        "xiHR": sp.simplify(beta / (1 + beta)),
        "Z1": sp.simplify(1 + sp.Float("0.4") * beta**2),
    }

    # reverse: partial local inversion witnesses from coefficient side
    xiHR = sp.Float(coeff["xiHR"])
    Mpl2 = sp.Float(coeff["Mpl2"])
    muH2 = sp.Float(coeff["muH2"])
    lambdaH = sp.Float(coeff["lambdaH"])

    beta_rec = sp.simplify(xiHR / (1 - xiHR))
    A_rec = sp.simplify(Mpl2 / (1 + beta_rec))
    omega_rec = sp.sqrt(sp.simplify(muH2 / A_rec))
    eta_rec = sp.sqrt(sp.simplify(lambdaH * (1 + beta_rec) - 1))

    # differences (numeric)
    diff_map = {
        "beta_abs": float(abs(beta_rec - beta)),
        "A_abs": float(abs(A_rec - A)),
        "omega_abs": float(abs(sp.N(omega_rec) - omega)),
        "eta_abs": float(abs(sp.N(eta_rec) - eta)),
    }

    payload = {
        "checkpoint": "P1694_S644",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "full_chain": "K_strict -> coefficients -> full L_total (nonskeleton anchor) -> sector EOM obligations",
        "kernel": "K_strict(d)=cos(omega*d+phi)/(1+beta*d^eta)",
        "kernel_input": kin,
        "forward_coefficient_map_symbolic": {k: str(v) for k, v in fwd.items()},
        "full_lagrangian_density_anchor": p1691.get("full_lagrangian_density", {}),
        "reverse_local_inversion": {
            "beta_rec": str(beta_rec),
            "A_rec": str(A_rec),
            "omega_rec": str(sp.N(omega_rec)),
            "eta_rec": str(sp.N(eta_rec)),
            "abs_differences": diff_map,
            "local_pass": all(v < 1e-12 for v in diff_map.values()),
            "scope": "local_parameter_recovery_only",
        },
        "bidirectional_status": {
            "forward_export": "EXPORTED",
            "reverse_export": "PARTIAL_LOCAL_ONLY",
            "global_variational_qg_theorem_closure": "KEEP_OPEN",
        },
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "open_obligations": [
            "full_covariant_EL_export_all_fields_nonproxy",
            "counterterm_flow_closure_theorem",
            "BRST_nilpotency_and_cohomology_theorem",
            "Cutkosky_full_sector_unitarity_theorem",
            "background_independence_family_theorem",
        ],
        "next_honest_step": "Podnieść reverse z local_parameter_recovery do pełnego odwrócenia: EOM-family -> integrable variational origin dla pełnego L_total + theorem-level QG closures.",
        "lay_summary": "Pokazaliśmy, że z kernela strict da się dojść do współczynników i pełnego lagranżianu oraz częściowo wrócić z parametrów. Pełny dowód całej kwantowej grawitacji nadal pozostaje otwarty.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
