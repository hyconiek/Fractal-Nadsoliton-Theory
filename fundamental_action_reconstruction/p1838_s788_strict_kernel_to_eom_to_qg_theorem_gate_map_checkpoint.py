#!/usr/bin/env python3
"""P1838 S788 strict kernel->EOM->QG theorem gate map checkpoint."""

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
    p1836 = load("p1836_s786_strict_full_lagrangian_non_skeleton_manifest_checkpoint.json")
    p1835 = load("p1835_s785_strict_forward_reverse_unlock_matrix_checkpoint.json")

    sectors = p1836.get("full_lagrangian_non_skeleton_manifest", {}).get("L_total", {})
    reverse_unlock = p1835.get("reverse_unlock_matrix", [])

    theorem_gate_map = {
        "TG1_BW": {
            "requires_forward": [
                "gravity_sector",
                "gauge_sector",
                "higgs_sector",
                "covariant_structures",
            ],
            "requires_reverse": ["helmholtz_integrability_global"],
            "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        },
        "TG2_BRST": {
            "requires_forward": [
                "gauge_sector",
                "fermion_sector",
                "interaction_terms",
                "covariant_structures",
            ],
            "requires_reverse": ["brst_nilpotency"],
            "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        },
        "TG3_CUT": {
            "requires_forward": [
                "fermion_sector",
                "interaction_terms",
                "curvature_corrections",
            ],
            "requires_reverse": ["cutkosky_unitarity", "renormalization_counterterm_closure", "background_independence"],
            "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        },
    }

    out = {
        "packet_id": "P1838",
        "stage_id": "S788",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "kernel_anchor": p1836.get("kernel_anchor", "K_strict_anchor_missing"),
        "full_sector_keys": list(sectors.keys()),
        "reverse_unlock_matrix": reverse_unlock,
        "theorem_gate_map": theorem_gate_map,
        "technical_progress": "Kernel->full Lagrangian sectors and reverse theorem obligations are now mapped directly to TG1/TG2/TG3 gate requirements.",
        "proven": "Each theorem gate requires both forward non-skeleton sectors and reverse theorem witnesses; neither side alone is sufficient.",
        "open": "All TG gates remain open pending missing sector exports and reverse-QG witness theorems.",
        "false_pass_risk": "Promoting TG gates from only forward or only reverse partials would violate strict bidirectional closure logic.",
        "next_honest_step": "Complete sector export gates then execute reverse theorem witnesses in unlock order and re-evaluate TG map.",
        "lay_explanation": "To mapa warunków bramek teorii: każda bramka potrzebuje i pełnych składników Lagrangianu, i dowodów po stronie QG.",
    }

    path = GEN / "p1838_s788_strict_kernel_to_eom_to_qg_theorem_gate_map_checkpoint.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
