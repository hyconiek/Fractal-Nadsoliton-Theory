#!/usr/bin/env python3
"""P1861 S811 strict B1 physical-pole projection and discontinuity certificate checkpoint."""

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
    p1860 = load("p1860_s810_strict_b1_cutkosky_kernel_sample_and_pole_residue_table_checkpoint.json")

    poles = ((p1860.get("pole_residue_seed_table") or {}).get("poles") or [])
    physical = [p for p in poles if float(p.get("s_pole", -1.0)) >= 0.0 and bool(p.get("nonnegative", False))]

    projection = {
        "rule": "retain poles with s_pole>=0 and residue>=0; discard ghost/negative-residue poles",
        "input_pole_count": len(poles),
        "physical_pole_count": len(physical),
        "physical_poles": physical,
    }

    # Seed certificate on projected subset (still not full theorem).
    disc_certificate = {
        "channel": "graviton->gauge_gauge",
        "statement": "Disc M_phys(s) = sum_{physical poles} Integral dPi_gg K_cut_phys >= 0 (seed projected corridor)",
        "projected_corridor": "s in [0.5, 8.0] from current seed grid",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "missing": [
            "full dressed propagator pole computation",
            "exact K_cut_phys phase-space integral",
            "background-family continuation beyond seed corridor",
        ],
    }

    out = {
        "packet_id": "P1861",
        "stage_id": "S811",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {
            "p1860_present": "pole_residue_seed_table" in p1860,
        },
        "physical_pole_projection": projection,
        "discontinuity_certificate_seed": disc_certificate,
        "strict_chain_extension": "K_strict -> K_cut sample -> pole table -> physical-pole projection -> projected discontinuity certificate",
        "proven": "Physical-pole projection rule and projected discontinuity certificate seed are explicit in strict B1 lane.",
        "open": "Projected certificate remains below theorem level until exact poles and integrals are computed.",
        "false_pass_risk": "Seed projection cannot be treated as final unitarity or TG3 closure proof.",
        "next_honest_step": "Replace seed poles with computed dressed poles and evaluate exact projected discontinuity integral with uncertainty-controlled positivity report.",
        "lay_explanation": "Odfiltrowaliśmy nie-fizyczne bieguny i zapisaliśmy test unitarności dla fizycznych, ale to nadal etap przed pełnym ścisłym dowodem.",
    }

    path = GEN / "p1861_s811_strict_b1_physical_pole_projection_and_disc_certificate_checkpoint.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
