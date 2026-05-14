from __future__ import annotations

import json
from pathlib import Path


def o_split_strict_contract(state_id: str) -> dict[str, object]:
    """Minimal strict-only O_split^(strict) I/O contract.

    Returns exactly three output objects and intentionally keeps selector closure
    open under QW-2191 unless a new strict-core selector source is provided.
    """
    return {
        "L_SM_channel": {
            "state_id": state_id,
            "mode": "strict_operational_projection",
        },
        "L_GR_carrier": {
            "state_id": state_id,
            "mode": "strict_geometric_carrier_projection",
        },
        "selector_status": "selector_source_missing",
    }


def main() -> None:
    out_dir = Path(__file__).resolve().parent / "generated"
    out_dir.mkdir(parents=True, exist_ok=True)

    route_state = "F_nadsoliton_to_LSM_plus_LGR"
    contract_output = o_split_strict_contract(state_id=route_state)

    summary = {
        "checkpoint": "P1521_S471",
        "status": "PASS_MINIMAL_STRICT_OSPLIT_IO_CONTRACT",
        "date_utc": "2026-05-13",
        "route": route_state,
        "strict_only": True,
        "legacy_bridge_used": False,
        "global_closure_claimed": False,
        "qw2191_closed": False,
        "contract_output": contract_output,
        "missing_strict_core_objects": [
            "strict_internal_selector_source_object",
            "selector_symmetry_breaking_premise_or_theorem",
            "strict_selector_provenance_witness",
        ],
    }

    out_path = out_dir / "p1521_s471_minimal_strict_osplit_io_contract_summary.json"
    out_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    print(f"[P1521] wrote {out_path}")


if __name__ == "__main__":
    main()
