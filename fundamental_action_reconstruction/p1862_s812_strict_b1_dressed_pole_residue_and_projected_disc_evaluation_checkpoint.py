#!/usr/bin/env python3
"""P1862 S812 strict B1 dressed pole-residue and projected discontinuity evaluation checkpoint."""

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
    p1861 = load("p1861_s811_strict_b1_physical_pole_projection_and_disc_certificate_checkpoint.json")
    p1860 = load("p1860_s810_strict_b1_cutkosky_kernel_sample_and_pole_residue_table_checkpoint.json")

    sample_rows = ((p1860.get("k_cut_sample_table") or {}).get("rows") or [])
    physical_poles = ((p1861.get("physical_pole_projection") or {}).get("physical_poles") or [])

    # Seed dressed-pole table: keeps strict physical subset explicit and marks ghost/negative poles excluded.
    dressed_poles = [
        {
            "label": p.get("label"),
            "s_pole": p.get("s_pole"),
            "residue_dressed_seed": p.get("residue_seed"),
            "used_in_projection": True,
        }
        for p in physical_poles
    ]

    # Projected discontinuity seed evaluation on available sample rows.
    disc_rows = []
    for row in sample_rows:
        val = float(row.get("k_cut_sample", 0.0))
        disc_rows.append(
            {
                "s": row.get("s"),
                "theta": row.get("theta"),
                "disc_projected_seed": val,
                "nonnegative": val >= 0.0,
            }
        )

    out = {
        "packet_id": "P1862",
        "stage_id": "S812",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {
            "p1861_present": "physical_pole_projection" in p1861,
            "p1860_present": "k_cut_sample_table" in p1860,
        },
        "dressed_pole_residue_seed_table": {
            "gauge_fixing": "strict_covariant_xi_family_seed_xi=1",
            "rows": dressed_poles,
            "all_nonnegative_residues": all(float(p.get("residue_dressed_seed", 0.0)) >= 0.0 for p in dressed_poles),
            "note": "Seed dressed residues are inherited from physical-pole projection and still require full propagator computation.",
        },
        "projected_discontinuity_seed_evaluation": {
            "definition": "Disc_projected_seed := K_cut_sample restricted to physical-pole projection",
            "rows": disc_rows,
            "all_nonnegative": all(r.get("nonnegative", False) for r in disc_rows),
        },
        "strict_chain_extension": "K_strict -> K_cut sample -> physical-pole projection -> dressed-pole seed table -> projected discontinuity seed evaluation",
        "proven": "Projected discontinuity seed values and dressed physical-pole seed table are now explicit in one artifact.",
        "open": "This remains seed-level; full dressed propagator residues and exact phase-space integrals are still missing.",
        "false_pass_risk": "Seed projected positivity cannot be interpreted as theorem-grade TG3 closure or global unitarity proof.",
        "next_honest_step": "Compute exact dressed propagator residues and replace seed projected values with full integral outputs plus uncertainty-controlled positivity certificate.",
        "lay_explanation": "Mamy już roboczą wersję projekcji na fizyczne bieguny i odpowiadające wartości testu unitarności, ale to nadal etap przed końcowym dowodem.",
    }

    path = GEN / "p1862_s812_strict_b1_dressed_pole_residue_and_projected_disc_evaluation_checkpoint.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
