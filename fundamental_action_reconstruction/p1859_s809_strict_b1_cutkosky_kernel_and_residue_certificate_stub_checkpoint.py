#!/usr/bin/env python3
"""P1859 S809 strict B1 Cutkosky kernel and residue certificate stub checkpoint."""

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
    p1858 = load("p1858_s808_strict_b1_cutkosky_seed_grid_evaluation_checkpoint.json")
    p1857 = load("p1857_s807_strict_b1_triangle_amplitude_seed_and_k5_instantiation_checkpoint.json")

    kernel_contract = {
        "channel": "graviton->gauge_gauge",
        "integral_kernel": "Disc M(s) = Integral dPi_gg K_cut(s,theta,phi; c_gr_i^ren, y_f, g_J)",
        "phase_space_domain": "2-body gauge phase-space in d=4-2*epsilon with MSbar_B1_seed",
        "matching_rule": "K_cut must reduce to seed proxy sign in low-order limit used by P1858",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    residue_certificate = {
        "target": "physical poles: Res D_phys(s_pole) >= 0 and ghost poles excluded by strict gauge-fixing constraints",
        "required_exports": [
            "dressed_graviton_propagator_pole_list",
            "residue_values_per_pole",
            "ghost_sector_exclusion_trace",
        ],
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    out = {
        "packet_id": "P1859",
        "stage_id": "S809",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {
            "p1858_present": "cutkosky_seed_grid" in p1858,
            "p1857_present": "cutkosky_discontinuity_numeric_symbolic_contract" in p1857,
        },
        "cutkosky_kernel_contract": kernel_contract,
        "residue_certificate_stub": residue_certificate,
        "strict_chain_extension": "K_strict -> renormalized coefficients -> explicit Cutkosky kernel -> residue positivity certificate",
        "proven": "Explicit discontinuity kernel contract and residue-certificate target are now formalized in strict B1 lane.",
        "open": "Kernel evaluation and residue data export remain missing.",
        "false_pass_risk": "Without explicit kernel integrals and residue tables no TG3 unitarity closure can be claimed.",
        "next_honest_step": "Export first computed K_cut sample table and corresponding dressed propagator pole/residue artifact under one fixed gauge-fixing choice.",
        "lay_explanation": "Mamy już zapis jak dokładnie ma wyglądać pełna całka unitarności i co trzeba udowodnić o biegunach, ale wciąż brakuje policzonych danych.",
    }

    path = GEN / "p1859_s809_strict_b1_cutkosky_kernel_and_residue_certificate_stub_checkpoint.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
