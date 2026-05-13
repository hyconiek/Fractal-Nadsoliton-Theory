#!/usr/bin/env python3
"""P1511 S4.61: admissible-perturbation boundary map for strict F=>LSM+LGR."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

P1500 = GEN / "p1500_s450_qw2191_selector_source_and_f_mapping_witness_summary.json"
SUMMARY = GEN / "p1511_s461_admissible_perturbation_boundary_map_f_to_lsm_lgr_summary.json"


def assumptions(selector: float, strict_internal: bool, w_sm: float, w_gr: float, o1: str, o2: str) -> dict[str, bool]:
    return {
        "A1_selector_positive": selector > 0.0,
        "A2_selector_strict_internal": strict_internal,
        "A3_weight_normalization": abs((w_sm + w_gr) - 1.0) <= 1e-9,
        "A4_shared_orientation": bool(o1) and o1 == o2,
        "A5_no_legacy_bridge": True,
    }


def main() -> None:
    s1500 = json.loads(P1500.read_text(encoding="utf-8"))
    s_internal = s1500.get("objects", {}).get("S_internal_v1", {})
    fmap = s1500.get("objects", {}).get("W_Fmap_v1", {})

    base_selector = float(s_internal.get("value", 0.0))
    strict_internal = bool(s_internal.get("strict_internal", False))
    base_w_sm = float(fmap.get("F_to_LSM_weight", 0.0))
    base_w_gr = float(fmap.get("F_to_LGR_weight", 0.0))
    base_o = str(s_internal.get("orientation", ""))
    base_fo = str(fmap.get("selector_orientation", ""))

    selector_scales = [0.85, 0.9, 1.0, 1.1, 1.15]
    weight_shifts = [-0.04, -0.02, 0.0, 0.02, 0.04]

    robust = []
    rejection = []

    for s in selector_scales:
        for dw in weight_shifts:
            w_sm = base_w_sm + dw
            w_gr = base_w_gr - dw
            aid = assumptions(base_selector * s, strict_internal, w_sm, w_gr, base_o, base_fo)
            row = {
                "selector_scale": s,
                "delta_w_sm": dw,
                "w_sm": w_sm,
                "w_gr": w_gr,
                "assumptions": aid,
                "pass": all(aid.values()),
            }
            (robust if row["pass"] else rejection).append(row)

    # explicit orientation boundary probe
    flip = assumptions(base_selector, strict_internal, base_w_sm, base_w_gr, "GR_preferred", base_fo)
    rejection.append({
        "selector_scale": 1.0,
        "delta_w_sm": 0.0,
        "w_sm": base_w_sm,
        "w_gr": base_w_gr,
        "assumptions": flip,
        "pass": all(flip.values()),
        "boundary_tag": "orientation_flip_boundary",
    })

    summary = {
        "packet": "P1511",
        "status": "PASS_ADMISSIBLE_BOUNDARY_MAP_EXPORTED",
        "scope": "STRICT_ONLY_F_NADSOLITON_TO_LSM_PLUS_LGR",
        "grid": {
            "selector_scales": selector_scales,
            "weight_shifts": weight_shifts,
        },
        "counts": {
            "robust_points": len(robust),
            "rejection_points": len(rejection),
        },
        "robust_zone": robust,
        "rejection_zone": rejection,
        "qw2191_closed": False,
        "legacy_bridge_used": False,
        "next_honest_step": "Implement P1512 minimal-condition core extraction from robust-zone points for strict coupled theorem strengthening.",
        "layman_explanation": "Narysowaliśmy mapę: które drobne zmiany model wytrzymuje, a które go łamią. Dzięki temu wiemy, gdzie teoria jest stabilna, a gdzie trzeba ją jeszcze wzmocnić.",
    }

    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1511] status={summary['status']} robust={len(robust)} rejection={len(rejection)}")


if __name__ == "__main__":
    main()
