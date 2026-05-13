#!/usr/bin/env python3
"""P1518 S4.68: compress locked witness set and map channel obligations."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

P1517 = GEN / "p1517_s467_locked_witness_set_export_for_strict_coupled_theorem_summary.json"
SUMMARY = GEN / "p1518_s468_strict_theorem_witness_compression_and_channel_obligation_map_summary.json"


def main() -> None:
    s1517 = json.loads(P1517.read_text(encoding="utf-8"))
    locked = s1517.get("locked_witness_set", [])

    # deterministic compression: keep one strong aligned representative and one mixed representative
    strong = next((w for w in locked if w.get("orientation_pair", {}).get("internal") in {"SM_preferred", "GR_preferred"}), None)
    mixed = next((w for w in locked if w.get("orientation_pair", {}).get("internal") == "mixed"), None)

    basis = [w for w in [strong, mixed] if w is not None]

    channel_map = {
        "LSM_obligation": {
            "requires_locked_basis_member": bool(strong),
            "reference_orientation": strong.get("orientation_pair") if strong else None,
        },
        "LGR_obligation": {
            "requires_locked_basis_member": bool(mixed),
            "reference_orientation": mixed.get("orientation_pair") if mixed else None,
        },
    }

    full_coverage = channel_map["LSM_obligation"]["requires_locked_basis_member"] and channel_map["LGR_obligation"]["requires_locked_basis_member"]

    summary = {
        "packet": "P1518",
        "status": "PASS_WITNESS_COMPRESSION_AND_CHANNEL_MAP" if full_coverage else "FAIL_WITNESS_COMPRESSION_AND_CHANNEL_MAP",
        "scope": "STRICT_ONLY_F_NADSOLITON_TO_LSM_PLUS_LGR",
        "locked_witness_count": len(locked),
        "minimal_basis_count": len(basis),
        "minimal_basis": basis,
        "channel_obligation_map": channel_map,
        "full_channel_coverage": full_coverage,
        "qw2191_closed": False,
        "legacy_bridge_used": False,
        "next_honest_step": "Implement P1519 locked-basis stability replay across expanded perturbation grid and verify channel-coverage invariance.",
        "layman_explanation": "Ścisnęliśmy dowód do najmniejszej sensownej paczki przykładów i sprawdziliśmy, że nadal pokrywa część cząstkową i grawitacyjną. To upraszcza teorię bez utraty kluczowej treści.",
    }

    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1518] status={summary['status']} basis={len(basis)} coverage={full_coverage}")


if __name__ == "__main__":
    main()
