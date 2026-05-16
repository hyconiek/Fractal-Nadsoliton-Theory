#!/usr/bin/env python3
"""P1860 S810 strict B1 Cutkosky kernel sample and pole/residue table checkpoint."""

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
    p1859 = load("p1859_s809_strict_b1_cutkosky_kernel_and_residue_certificate_stub_checkpoint.json")
    p1853 = load("p1853_s803_strict_b1_symbolic_evaluation_and_scheme_trace_checkpoint.json")

    coeffs = ((p1853.get("b1_symbolic_evaluation") or {}).get("evaluated_coefficients") or {})
    a_r2 = float((coeffs.get("a_R2") or {}).get("numeric_20d", "0") or 0.0)
    a_ric2 = float((coeffs.get("a_Ric2") or {}).get("numeric_20d", "0") or 0.0)

    # Seed sample of explicit kernel values (still proxy-level, but explicit table export).
    s_points = [0.5, 1.0, 2.0, 4.0]
    angle_points = [0.0, 0.5, 1.0]
    k_cut_rows = []
    for s in s_points:
        for th in angle_points:
            k_val = (a_r2 + a_ric2 * (1.0 + th * th)) * s
            k_cut_rows.append({"s": s, "theta": th, "k_cut_sample": k_val, "nonnegative": k_val >= 0.0})

    pole_residue_table = {
        "gauge_fixing": "strict_covariant_xi_family_seed_xi=1",
        "poles": [
            {"label": "spin2_physical", "s_pole": 0.0, "residue_seed": 1.0, "nonnegative": True},
            {"label": "massive_seed_mode", "s_pole": 1.0, "residue_seed": 0.12, "nonnegative": True},
            {"label": "ghost_candidate", "s_pole": -1.0, "residue_seed": -0.05, "nonnegative": False},
        ],
        "ghost_exclusion_seed_rule": "only poles in physical s>=0 corridor with nonnegative residue can enter unitarity channel sum",
    }

    out = {
        "packet_id": "P1860",
        "stage_id": "S810",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {
            "p1859_present": "cutkosky_kernel_contract" in p1859,
            "p1853_present": "b1_symbolic_evaluation" in p1853,
        },
        "k_cut_sample_table": {
            "definition": "K_cut_sample(s,theta) := (a_R2 + a_Ric2*(1+theta^2))*s",
            "rows": k_cut_rows,
            "all_nonnegative": all(r["nonnegative"] for r in k_cut_rows),
        },
        "pole_residue_seed_table": pole_residue_table,
        "strict_chain_extension": "K_strict -> evaluated coefficients -> explicit K_cut sample table -> pole/residue seed classification",
        "proven": "First explicit K_cut sample table and pole/residue seed export are now available for strict B1 TG3 preparation.",
        "open": "Sample-level kernel and residue seeds are not sufficient for full discontinuity integral theorem closure.",
        "false_pass_risk": "Seed pole table with proxy residues cannot be interpreted as final ghost-free unitarity proof.",
        "next_honest_step": "Replace seed pole residues with computed dressed propagator residues and evaluate full phase-space discontinuity integral with certified physical-pole projection.",
        "lay_explanation": "Mamy już przykładową tabelę wartości jądra unitarności i biegunów, ale to nadal wersja robocza przed pełnym dowodem.",
    }

    path = GEN / "p1860_s810_strict_b1_cutkosky_kernel_sample_and_pole_residue_table_checkpoint.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
