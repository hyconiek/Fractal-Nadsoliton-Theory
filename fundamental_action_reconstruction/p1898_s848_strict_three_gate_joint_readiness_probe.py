#!/usr/bin/env python3
"""P1898 S848 strict three-gate joint readiness probe (G_R, G_U, G_BI)."""
from __future__ import annotations
import json
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def load(name: str) -> dict:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def main() -> None:
    p1894 = load("p1894_s844_strict_gr_promotion_readiness_update_probe.json")
    p1896 = load("p1896_s846_strict_gu_promotion_readiness_probe.json")
    p1897 = load("p1897_s847_strict_gbi_preclosure_transport_bridge_probe.json")

    beta, eta, omega, phi = sp.symbols("beta eta omega phi", positive=True, real=True)
    lam, y, kappa2 = sp.symbols("lam y kappa2", positive=True, real=True)
    eps, mu2, s = sp.symbols("eps mu2 s", positive=True, real=True)

    c0 = sp.cos(omega + phi) / (beta + 1)
    c1 = (-beta * eta * sp.cos(omega + phi) - omega * (beta + 1) * sp.sin(omega + phi)) / (beta + 1) ** 2
    lam_eff = lam * (1 + c1**2)
    y_eff = y * (1 + c0 / 2)

    I_lam = 1 / (16 * sp.pi**2) * (1 / eps + sp.log(mu2 / s))
    I_y = 1 / (16 * sp.pi**2) * (1 / eps + sp.log(mu2 / s) + sp.Rational(1, 2))

    G_R = lam_eff * I_lam
    G_U = y_eff**2 * I_y
    G_BI = (y_eff**2 + kappa2) * I_y

    r_gr, r_gu, r_gbi = sp.symbols("r_gr r_gu r_gbi", nonnegative=True, real=True)
    joint_score = sp.simplify((r_gr + r_gu + r_gbi) / 3)

    out = {
        "packet_id": "P1898",
        "stage_id": "S848",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {
            "p1894_present": "promotion_readiness" in p1894,
            "p1896_present": "promotion_readiness" in p1896,
            "p1897_present": "promotion_readiness" in p1897,
        },
        "strict_chain_step": "K_strict -> coefficients -> full non-skeleton L_SM+L_GR -> EOM -> three-gate joint readiness",
        "gate_signals": {"G_R": str(sp.expand(G_R)), "G_U": str(sp.expand(G_U)), "G_BI": str(sp.expand(G_BI))},
        "joint_readiness": {
            "score_symbol": str(joint_score),
            "components": {
                "r_gr": "renormalization-gate readiness",
                "r_gu": "unitarity-gate readiness",
                "r_gbi": "background-independence-gate readiness",
            },
            "strict_promotion_contract": "joint score >= 1 and each component backed by witness-grade evidence",
        },
        "strict_core_closure_missing_items": {
            "renormalization": "Need G_R witness package finalized.",
            "unitarity": "Need G_U discontinuity witness package finalized.",
            "background_independence": "Need G_BI transport theorem/witness finalized.",
            "selector_obstruction": "QW-2191 remains open.",
        },
        "false_pass_guard": "Joint readiness score is not strict-core closure proof.",
        "next_honest_step": "Promote first gate with full evidence, then recompute joint score and iterate until all three gates are witness-closed.",
        "lay_explanation": "To wspólny licznik postępu trzech kluczowych testów teorii: renormalizacji, unitarności i niezależności od tła.",
    }

    path = GEN / "p1898_s848_strict_three_gate_joint_readiness_probe.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
