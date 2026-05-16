#!/usr/bin/env python3
"""P1897 S847 strict G_BI preclosure transport bridge probe."""
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
    p1896 = load("p1896_s846_strict_gu_promotion_readiness_probe.json")
    p1887 = load("p1887_s837_strict_nontrivial_transport_branch_probe.json")

    beta, eta, omega, phi = sp.symbols("beta eta omega phi", positive=True, real=True)
    y, kappa2 = sp.symbols("y kappa2", positive=True, real=True)
    eps, mu2, s = sp.symbols("eps mu2 s", positive=True, real=True)
    sigma2 = sp.symbols("sigma2", nonnegative=True, real=True)

    c0 = sp.cos(omega + phi) / (beta + 1)
    y_eff = y * (1 + c0 / 2)
    I_y = 1 / (16 * sp.pi**2) * (1 / eps + sp.log(mu2 / s) + sp.Rational(1, 2))

    G_BI = (y_eff**2 + kappa2) * I_y

    nu = sp.symbols("nu", real=True)
    # Bridge from FRW to Bianchi-I for BI gate carrier.
    G_BI_frw = G_BI
    G_BI_b1 = G_BI * (1 + nu * sigma2)
    delta_bi = sp.expand(G_BI_b1 - G_BI_frw)

    b_cov, b_overlap, b_scheme = sp.symbols("b_cov b_overlap b_scheme", nonnegative=True, real=True)
    readiness_score = sp.simplify((b_cov + b_overlap + b_scheme) / 3)

    out = {
        "packet_id": "P1897",
        "stage_id": "S847",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {
            "p1896_present": "promotion_readiness" in p1896,
            "p1887_present": "nontrivial_transport_branch" in p1887,
        },
        "strict_chain_step": "K_strict -> coefficients -> full L_total transport carriers -> G_BI preclosure bridge",
        "effective_coefficients": {"y_eff": str(y_eff)},
        "gbi_bridge_quantities": {
            "G_BI_frw": str(sp.expand(G_BI_frw)),
            "G_BI_b1": str(sp.expand(G_BI_b1)),
            "delta_bi": str(delta_bi),
            "branch_parameter": "nu",
        },
        "promotion_readiness": {
            "score_symbol": str(readiness_score),
            "score_components": {
                "b_cov": "covariant transport equation coverage",
                "b_overlap": "FRW/Bianchi atlas overlap coherence",
                "b_scheme": "same renormalization-scheme continuity across backgrounds",
            },
            "promotion_threshold_contract": "score_symbol >= 1 with theorem-grade transport witness",
        },
        "strict_core_closure_missing_items": {
            "background_independence": "Need nontrivial nu-branch witness proving controlled delta_bi behavior.",
            "unitarity": "Need compatibility with G_U discontinuity closures.",
            "renormalization": "Need shared finite-part lock across FRW and Bianchi-I branches.",
            "selector_obstruction": "QW-2191 remains open.",
        },
        "false_pass_guard": "G_BI bridge probe is preclosure and not a background-independence theorem.",
        "next_honest_step": "Attach explicit transport witness data for nu!=0 branch and attempt first controlled G_BI readiness upgrade.",
        "lay_explanation": "Po renormalizacji i unitarności przychodzi krok tła: sprawdzamy, czy teoria działa spójnie także między różnymi geometriami czasoprzestrzeni.",
    }

    path = GEN / "p1897_s847_strict_gbi_preclosure_transport_bridge_probe.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
