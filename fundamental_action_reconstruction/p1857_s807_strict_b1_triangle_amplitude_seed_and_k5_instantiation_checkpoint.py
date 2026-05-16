#!/usr/bin/env python3
"""P1857 S807 strict B1 triangle amplitude seed and k5 instantiation checkpoint."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def load(name: str) -> dict:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def main() -> None:
    p1856 = load("p1856_s806_strict_b1_representation_and_k5_triangle_seed_checkpoint.json")
    p1853 = load("p1853_s803_strict_b1_symbolic_evaluation_and_scheme_trace_checkpoint.json")

    rep_seed = (p1856.get("strict_representation_seed") or {}).get("fermion_representations_seed", [])
    a_ric2 = ((p1853.get("b1_symbolic_evaluation") or {}).get("evaluated_coefficients") or {}).get("a_Ric2", {}).get("symbolic", "a_Ric2_B1")

    triangle_seed = {
        "scheme": "MSbar_B1_seed",
        "amplitude_form": "I_triangle_f = (1/(16*pi^2))*C_triangle_f*log(mu_RG^2/m_f^2) + finite_f",
        "strict_coupling_link": f"C_triangle_f := strict_triangle_norm_f * ({a_ric2})",
        "regularization": "dimensional_regularization_d=4-2*epsilon",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    k5_instantiation = {
        "definition": "k5 = sum_f Tr[T_a {T_b,T_c}] * Y_f * I_triangle_f",
        "family_sum_seed_count": len(rep_seed),
        "seed_partial_sum_note": "Current artifact exports first-family seed entries; full 3-family/chiral-complete sum required.",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    discontinuity_numeric_symbolic_contract = {
        "channel": "graviton->gauge_gauge",
        "disc_form": "Disc M(s)=Integral dPi_gg M(grav->gg)M*(grav->gg)",
        "numeric_symbolic_target": "evaluate Disc M(s) at seed grid s in B1 and verify nonnegative sign",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    out = {
        "packet_id": "P1857",
        "stage_id": "S807",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {
            "p1856_present": "k5_triangle_seed" in p1856,
            "p1853_present": "b1_symbolic_evaluation" in p1853,
        },
        "triangle_amplitude_seed": triangle_seed,
        "k5_instantiation_seed": k5_instantiation,
        "cutkosky_discontinuity_numeric_symbolic_contract": discontinuity_numeric_symbolic_contract,
        "proven": "Regularized triangle amplitude seed and k5 instantiation formula are now explicit in B1 strict scheme.",
        "open": "Full family-complete k5 value and explicit discontinuity evaluation/positivity certificate remain open.",
        "false_pass_risk": "Seed instantiation without complete sums and evaluated integrals does not discharge BRST/Cutkosky gates.",
        "next_honest_step": "Complete 3-family chiral sum for k5 and run first seed-grid Disc M(s) evaluation with signed positivity report.",
        "lay_explanation": "Mamy już roboczy wzór na brakującą amplitudę trójkątną i sposób liczenia k5, ale jeszcze trzeba policzyć pełną sumę i sprawdzić znak całki unitarności.",
    }

    path = GEN / "p1857_s807_strict_b1_triangle_amplitude_seed_and_k5_instantiation_checkpoint.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
