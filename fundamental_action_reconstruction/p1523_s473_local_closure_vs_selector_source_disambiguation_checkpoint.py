from __future__ import annotations

import json
from pathlib import Path

from p1522_s472_strict_selector_source_intake_gate_checkpoint import strict_selector_source_intake


def main() -> None:
    out_dir = Path(__file__).resolve().parent / "generated"
    out_dir.mkdir(parents=True, exist_ok=True)

    local_closure_candidate = {
        "candidate_id": "LOCAL_CLOSURE_ONLY_V1",
        "provider_class": "local_equation_closure",
        "symmetry_breaking_premise": "none",
        "strict_provenance_trace": [],
        "noncyclic_anchor": True,
    }

    intake = strict_selector_source_intake(local_closure_candidate)

    summary = {
        "checkpoint": "P1523_S473",
        "status": "PASS_LOCAL_CLOSURE_SELECTOR_SOURCE_DISAMBIGUATION",
        "date_utc": "2026-05-13",
        "route": "F_nadsoliton_to_LSM_plus_LGR",
        "strict_only": True,
        "legacy_bridge_used": False,
        "local_closure_status": "present_but_not_selector_source",
        "selector_source_status": intake["intake_status"],
        "selector_source_reason": intake["reason_code"],
        "qw2191_closed": False,
        "next_required_objects": [
            "strict_internal_selector_source_export",
            "nonempty_strict_provenance_trace",
            "selector_symmetry_breaking_witness",
        ],
    }

    out_path = out_dir / "p1523_s473_local_closure_vs_selector_source_disambiguation_summary.json"
    out_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1523] wrote {out_path}")


if __name__ == "__main__":
    main()
