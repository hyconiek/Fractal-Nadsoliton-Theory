#!/usr/bin/env python3
"""P1870 S820 strict FRW-background metric residual probe (componentwise fill step)."""
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
    p1869 = load("p1869_s819_strict_flat_background_component_residual_probe.json")

    t = sp.symbols("t", real=True)
    H, kappa2, Lambda_cc, rho, p = sp.symbols("H kappa2 Lambda_cc rho p", real=True)

    # FRW with constant H (de Sitter-like test slice): component Einstein tensor entries.
    G00 = 3 * H**2
    Gii = -3 * H**2

    # Perfect-fluid diagonal stress tensor in comoving frame.
    T00 = rho
    Tii = p

    Rg00 = sp.simplify(G00 + Lambda_cc - kappa2 * T00)
    Rg11 = sp.simplify(Gii + Lambda_cc - kappa2 * Tii)

    # Isotropy gives 22,33 same class as 11 (up to scale factors, tracked as same reduced class here).
    out = {
        "packet_id": "P1870",
        "stage_id": "S820",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {"p1869_present": "residual_fill_result" in p1869},
        "declared_background_family": "frw_constant_H_slice",
        "assumptions": [
            "metric family: spatially flat FRW",
            "H(t)=const on probe slice",
            "T_mn = diag(rho, p, p, p)",
            "scalar/gauge local matter-source residual channels deferred to dedicated sector fills",
        ],
        "einstein_component_residual_probe": {
            "R_g_00": str(Rg00),
            "R_g_11_class": str(Rg11),
            "R_g_22_class": str(Rg11),
            "R_g_33_class": str(Rg11),
            "closure_conditions": {
                "friedmann_like": "3*H**2 + Lambda_cc - kappa2*rho = 0",
                "pressure_like": "-3*H**2 + Lambda_cc - kappa2*p = 0",
            },
            "component_status": "OPEN_OBSTRUCTION_WITH_TRACE",
            "obstruction_reason": "Residuals vanish only on constrained manifold linking (H, Lambda_cc, rho, p); generic FRW slice is not residual-zero.",
        },
        "bidirectional_note": {
            "forward": "Given strict kernel-derived couplings, closure requires background+matter satisfying Einstein residual constraints.",
            "reverse": "Any accepted FRW residual-zero branch imposes constraints to propagate back into admissible strict closure sector data.",
        },
        "qg_closure_gap": {
            "renormalization": "Need one-loop counterterm closure on same FRW-reduced operator basis.",
            "unitarity": "Need Cutkosky/optical-theorem witness compatible with this background family.",
            "background_independence": "Need atlas lift proving consistency beyond chosen FRW slice.",
        },
        "false_pass_guard": "FRW algebraic component probe is not theorem-grade global closure.",
        "next_honest_step": "Attach explicit strict coefficient sector to (rho,p) constitutive map and compute whether FRW closure manifold is nonempty under same renormalization scheme.",
        "lay_explanation": "Na tle kosmologicznym FRW dostajemy konkretne równania składowych grawitacji. One nie znikają same z siebie — trzeba spełnić warunki między ekspansją Wszechświata, energią i ciśnieniem.",
    }

    path = GEN / "p1870_s820_strict_frw_background_metric_residual_probe.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
