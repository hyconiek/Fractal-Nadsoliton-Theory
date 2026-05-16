#!/usr/bin/env python3
"""P1895 S845 strict G_U preclosure Cutkosky bridge probe."""
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
    p1876 = load("p1876_s826_strict_channel_amplitude_discontinuity_table_probe.json")

    beta, eta, omega, phi = sp.symbols("beta eta omega phi", positive=True, real=True)
    y, kappa2 = sp.symbols("y kappa2", positive=True, real=True)
    eps, mu2, s = sp.symbols("eps mu2 s", positive=True, real=True)

    c0 = sp.cos(omega + phi) / (beta + 1)
    y_eff = y * (1 + c0 / 2)

    I_y = 1 / (16 * sp.pi**2) * (1 / eps + sp.log(mu2 / s) + sp.Rational(1, 2))
    GU_ss = y_eff**2 * I_y
    GU_ff = (y_eff**2 + kappa2) * I_y

    rho2_ss, rho2_ff = sp.symbols("rho2_ss rho2_ff", nonnegative=True, real=True)
    ImM_ss, ImM_ff = sp.symbols("ImM_ss ImM_ff", real=True)

    Mss = y_eff**2 * (1 + y_eff**2 * I_y)
    Mff = (y_eff**2 + kappa2) * (1 + y_eff**2 * I_y)

    defect_ss = sp.expand(ImM_ss - rho2_ss * Mss**2)
    defect_ff = sp.expand(ImM_ff - rho2_ff * Mff**2)

    out = {
        "packet_id": "P1895",
        "stage_id": "S845",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {
            "p1894_present": "promotion_readiness" in p1894,
            "p1876_present": "cutkosky_discontinuity_table" in p1876,
        },
        "strict_chain_step": "K_strict -> coefficients -> full L_total Yukawa/fermion slots -> G_U bridge defects",
        "effective_coefficients": {"y_eff": str(y_eff)},
        "gu_bridge_quantities": {
            "GU_ss": str(sp.expand(GU_ss)),
            "GU_ff": str(sp.expand(GU_ff)),
            "Mss": str(sp.expand(Mss)),
            "Mff": str(sp.expand(Mff)),
        },
        "cutkosky_bridge_defects": {
            "defect_ss": str(defect_ss),
            "defect_ff": str(defect_ff),
            "target": "defect_ss=0 and defect_ff=0 under one shared renormalization scheme",
        },
        "strict_core_closure_missing_items": {
            "unitarity": "Need explicit ImM_ss/ImM_ff from diagrammatic discontinuities.",
            "renormalization": "Need finite-part lock compatibility with promoted G_R lane.",
            "background_independence": "Need branch transport of GU bridge on FRW/BianchiI.",
            "selector_obstruction": "QW-2191 remains open.",
        },
        "false_pass_guard": "G_U bridge probe is preclosure and not unitarity theorem closure.",
        "next_honest_step": "Attach explicit discontinuity integrals for ImM_ss/ImM_ff and attempt first G_U promotion readiness metric.",
        "lay_explanation": "Po renormalizacji przechodzimy do unitarności: sprawdzamy, czy amplitudy rozpraszania spełniają równania Cutkosky'ego w tym samym schemacie.",
    }

    path = GEN / "p1895_s845_strict_gu_preclosure_cutkosky_bridge_probe.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
