#!/usr/bin/env python3
"""P1896 S846 strict G_U promotion readiness probe."""
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
    p1895 = load("p1895_s845_strict_gu_preclosure_cutkosky_bridge_probe.json")
    p1894 = load("p1894_s844_strict_gr_promotion_readiness_update_probe.json")

    beta, eta, omega, phi = sp.symbols("beta eta omega phi", positive=True, real=True)
    y, kappa2 = sp.symbols("y kappa2", positive=True, real=True)
    eps, mu2, s = sp.symbols("eps mu2 s", positive=True, real=True)

    c0 = sp.cos(omega + phi) / (beta + 1)
    y_eff = y * (1 + c0 / 2)
    I_y = 1 / (16 * sp.pi**2) * (1 / eps + sp.log(mu2 / s) + sp.Rational(1, 2))

    GU_ss = y_eff**2 * I_y
    GU_ff = (y_eff**2 + kappa2) * I_y

    u_disc, u_scheme, u_transport = sp.symbols("u_disc u_scheme u_transport", nonnegative=True, real=True)
    readiness_score = sp.simplify((u_disc + u_scheme + u_transport) / 3)

    out = {
        "packet_id": "P1896",
        "stage_id": "S846",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {
            "p1895_present": "cutkosky_bridge_defects" in p1895,
            "p1894_present": "promotion_readiness" in p1894,
        },
        "strict_chain_step": "K_strict -> coefficients -> full L_total Yukawa/fermion sectors -> G_U promotion readiness",
        "effective_coefficients": {"y_eff": str(y_eff)},
        "gu_signal": {
            "GU_ss": str(sp.expand(GU_ss)),
            "GU_ff": str(sp.expand(GU_ff)),
        },
        "promotion_readiness": {
            "score_symbol": str(readiness_score),
            "score_components": {
                "u_disc": "explicit discontinuity (ImM) coverage",
                "u_scheme": "shared renormalization-scheme coherence with G_R lane",
                "u_transport": "FRW/BianchiI transport compatibility for unitarity channel",
            },
            "promotion_threshold_contract": "score_symbol >= 1 with witness-grade evidence",
        },
        "strict_core_closure_missing_items": {
            "unitarity": "Need explicit ImM_ss/ImM_ff discontinuities to close defects.",
            "renormalization": "Need consistency with promoted/near-promoted G_R dataset.",
            "background_independence": "Need transport witness for GU channels.",
            "selector_obstruction": "QW-2191 remains open.",
        },
        "false_pass_guard": "Readiness score does not imply G_U closure without discontinuity witnesses.",
        "next_honest_step": "Attach computed ImM discontinuities and attempt first controlled G_U status update on shared scheme.",
        "lay_explanation": "To odpowiednik kroku dla renormalizacji, ale dla unitarności: mierzymy, czy mamy już dość dowodów, by podnieść status bramy G_U.",
    }

    path = GEN / "p1896_s846_strict_gu_promotion_readiness_probe.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
