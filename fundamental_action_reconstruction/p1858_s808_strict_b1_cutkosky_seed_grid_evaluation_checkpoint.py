#!/usr/bin/env python3
"""P1858 S808 strict B1 Cutkosky seed-grid evaluation checkpoint."""

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
    p1857 = load("p1857_s807_strict_b1_triangle_amplitude_seed_and_k5_instantiation_checkpoint.json")
    p1853 = load("p1853_s803_strict_b1_symbolic_evaluation_and_scheme_trace_checkpoint.json")

    coeffs = ((p1853.get("b1_symbolic_evaluation") or {}).get("evaluated_coefficients") or {})
    a_r2 = float((coeffs.get("a_R2") or {}).get("numeric_20d", "0") or 0.0)
    a_ric2 = float((coeffs.get("a_Ric2") or {}).get("numeric_20d", "0") or 0.0)

    # Seed-only proxy (not theorem): positive corridor indicator over sample s-grid.
    s_grid = [0.5, 1.0, 2.0, 4.0, 8.0]
    rows = []
    for s in s_grid:
        disc_proxy = (a_r2 + 0.5 * a_ric2) * s
        rows.append({
            "s": s,
            "disc_proxy": disc_proxy,
            "nonnegative": disc_proxy >= 0.0,
        })

    all_nonnegative = all(r["nonnegative"] for r in rows)

    out = {
        "packet_id": "P1858",
        "stage_id": "S808",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {
            "p1857_present": "cutkosky_discontinuity_numeric_symbolic_contract" in p1857,
            "p1853_present": "b1_symbolic_evaluation" in p1853,
        },
        "cutkosky_seed_grid": {
            "channel": "graviton->gauge_gauge",
            "scheme": "MSbar_B1_seed",
            "proxy_definition": "Disc_proxy(s) := (a_R2 + 0.5*a_Ric2)*s",
            "grid_results": rows,
            "all_nonnegative_on_seed_grid": all_nonnegative,
        },
        "strict_chain_extension": "K_strict -> B1 coefficient evaluation -> first seed-grid discontinuity proxy scan",
        "proven": "First numeric-symbolic seed-grid scan for discontinuity sign has been exported in strict B1 lane.",
        "open": "Proxy positivity on seed grid is not a full Cutkosky integral theorem and does not prove global unitarity.",
        "false_pass_risk": "Seed-grid proxy scan cannot replace explicit discontinuity integral evaluation and residue analysis.",
        "next_honest_step": "Replace Disc_proxy with explicit discontinuity integral kernel and add residue/pole-structure certificate per channel.",
        "lay_explanation": "Sprawdziliśmy pierwszy uproszczony test znaku dla unitarności na kilku punktach, ale to jeszcze nie jest pełny dowód fizyczny.",
    }

    path = GEN / "p1858_s808_strict_b1_cutkosky_seed_grid_evaluation_checkpoint.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
