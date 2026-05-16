#!/usr/bin/env python3
"""P1856 S806 strict B1 representation and k5 triangle seed checkpoint."""

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
    p1855 = load("p1855_s805_strict_b1_gauge_fermion_k5_and_cutkosky_integral_stub_checkpoint.json")
    p1853 = load("p1853_s803_strict_b1_symbolic_evaluation_and_scheme_trace_checkpoint.json")

    reps = {
        "fermion_representations_seed": [
            {"family": "f1", "su3": "3", "su2": "2", "Y": "1/6"},
            {"family": "f1", "su3": "3", "su2": "1", "Y": "2/3"},
            {"family": "f1", "su3": "3", "su2": "1", "Y": "-1/3"},
            {"family": "f1", "su3": "1", "su2": "2", "Y": "-1/2"},
            {"family": "f1", "su3": "1", "su2": "1", "Y": "-1"},
        ],
        "note": "Seed table for strict B1 anomaly-triangle bookkeeping; full 3-family expansion required for theorem closure.",
    }

    eval_coeffs = ((p1853.get("b1_symbolic_evaluation") or {}).get("evaluated_coefficients") or {})
    a_ric2 = (eval_coeffs.get("a_Ric2") or {}).get("symbolic", "a_Ric2_B1")

    k5_triangle_seed = {
        "formula": "k5 = sum_f Tr[T_a {T_b,T_c}] * Y_f * I_triangle_f(B1)",
        "triangle_integral_seed": "I_triangle_f(B1) := C_triangle * a_Ric2_B1 (scheme-coupled seed)",
        "substitution_seed": {"a_Ric2_B1": a_ric2, "C_triangle": "strict_triangle_norm"},
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "trace": "Need explicit regulator-consistent triangle amplitude evaluation and full family sum with anomaly cancellation check.",
    }

    out = {
        "packet_id": "P1856",
        "stage_id": "S806",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {
            "p1855_present": "k5_gauge_fermion_contract" in p1855,
            "p1853_present": "b1_symbolic_evaluation" in p1853,
        },
        "strict_representation_seed": reps,
        "k5_triangle_seed": k5_triangle_seed,
        "strict_chain_extension": "K_strict -> B1 coeff eval -> representation seed -> k5 triangle seed -> BRST anomaly closure target",
        "proven": "Strict B1 representation seed and first k5 triangle coupling seed are now explicit.",
        "open": "Full anomaly-cancellation sum and regulator-stable triangle evaluation remain open.",
        "false_pass_risk": "Seed representation/triangle formulas are not a completed BRST anomaly cancellation theorem.",
        "next_honest_step": "Expand to all families/chiral sectors, compute regulator-stable triangle amplitudes, and export explicit A_B1 cancellation witness.",
        "lay_explanation": "Ustaliliśmy bazową tabelę cząstek i szkic wzoru na brakujący składnik anomalii, ale trzeba jeszcze policzyć pełny wynik dla wszystkich rodzin.",
    }

    path = GEN / "p1856_s806_strict_b1_representation_and_k5_triangle_seed_checkpoint.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
