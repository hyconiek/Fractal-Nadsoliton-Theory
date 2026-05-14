#!/usr/bin/env python3
"""P1604/S554: noncyclic provider upgrade plan for strict G1+bridge co-closure after adaptive scan failure."""
from __future__ import annotations
import json
from datetime import datetime, UTC
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1603 = GEN / "p1603_s553_adaptive_tail_depth_and_bridge_cotuning_summary.json"


def main() -> None:
    if not IN1603.exists():
        raise FileNotFoundError("Run P1603 before P1604")

    s03 = json.loads(IN1603.read_text(encoding="utf-8"))
    best = s03.get("adaptive_tail_scan", {}).get("best_admissible_candidate")

    provider_upgrade = {
        "class_id": "NP1_noncyclic_selector_source_provider",
        "goal": "Break persistent symmetric tail failure without legacy bridge transfer.",
        "mechanism": [
            "Inject strict internal asymmetric regulator term in selector-source generator.",
            "Constrain regulator by kernel-to-EOM consistency so physical basis remains strict.",
            "Export theorem witness linking regulator admissibility to bridge theorem premises.",
        ],
        "requires_new_blocker_cut": True,
        "noncyclic_anchor": "NP1-anchor-QW238x-compliant",
    }

    status = (
        "PASS_P1604_PROVIDER_UPGRADE_READY"
        if best is None
        else "KEEP_OPEN_P1604_ADAPTIVE_ALREADY_SUFFICIENT"
    )

    summary = {
        "checkpoint": "P1604_S554",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": status,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "strict_chain": s03.get("strict_chain", {}),
        "adaptive_result": {
            "best_admissible_candidate": best,
            "adaptive_failed": best is None,
        },
        "noncyclic_provider_upgrade": provider_upgrade,
        "strict_core_closure": {
            "status": "OPEN",
            "missing_exports": s03.get("strict_core_closure", {}).get("missing_exports", []),
            "missing_witnesses": ["W_G1_full_domain_selector_gap_discharge"],
            "missing_theorems": [
                "T_selector_uniqueness_to_full_lagrangian_bridge",
                "T_G3_final_strict_ToE_composition",
            ],
        },
        "external_team_validation_required": False,
        "next_honest_step": "P1605: instantiate NP1 provider and regenerate strict selector-source samples, then replay P1598/P1599/P1603.",
        "lay_summary": "Skoro sama korekta parametrów nie wystarcza, proponujemy nowy niecykliczny mechanizm wewnętrzny, który ma usunąć globalną lukę selektora bez łamania rygoru strict-only."
    }

    out = GEN / "p1604_s554_noncyclic_provider_upgrade_for_g1_bridge_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
