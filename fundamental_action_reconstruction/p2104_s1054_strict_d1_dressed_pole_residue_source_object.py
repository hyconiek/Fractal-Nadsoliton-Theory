#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2103 = GEN / "p2103_s1053_strict_dressed_discontinuity_backend_import_contract.json"
IN_1862 = GEN / "p1862_s812_strict_b1_dressed_pole_residue_and_projected_disc_evaluation_checkpoint.json"
OUT = GEN / "p2104_s1054_strict_d1_dressed_pole_residue_source_object.json"
MD = GEN / "p2104_s1054_strict_d1_dressed_pole_residue_source_object.md"

SCHEMA_VERSION = "p2104_s1054_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2103 = load(IN_2103)
    p1862 = load(IN_1862)

    pre_ready = p2103.get("result_kind") == "PASS_STRICT_DRESSED_DISCONTINUITY_BACKEND_IMPORT_CONTRACT_WITH_TRACE"

    dressed_seed = p1862.get("dressed_pole_residue_seed_table", {}) if isinstance(p1862, dict) else {}
    rows = dressed_seed.get("rows", []) if isinstance(dressed_seed, dict) else []

    has_rows = isinstance(rows, list) and len(rows) > 0

    result_kind = (
        "PASS_STRICT_D1_DRESSED_POLE_RESIDUE_SOURCE_OBJECT_WITH_TRACE__SEED_LEVEL"
        if pre_ready and has_rows
        else "OPEN_STRICT_D1_DRESSED_POLE_RESIDUE_SOURCE_OBJECT_BLOCKED"
    )

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2104",
        "stage_id": "S1054",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": result_kind,
        "depends_on": {
            "p2103_present": p2103.get("_missing") is None,
            "p1862_present": p1862.get("_missing") is None,
            "preconditions_ready": pre_ready,
            "seed_rows_present": has_rows,
        },
        "d1_source_object": {
            "object_name": "dressed_pole_residue_seed_source_on_common_basis",
            "channel": "graviton -> gauge_gauge",
            "rows_count": len(rows) if isinstance(rows, list) else 0,
            "source_artifact": "p1862_s812_strict_b1_dressed_pole_residue_and_projected_disc_evaluation_checkpoint.json",
            "scope_limit": "seed-level source object; not full dressed residue theorem object",
        },
        "backend_import_blocker_update": {
            "D1_dressed_pole_residue_backend_on_common_basis": "COMPUTED" if (pre_ready and has_rows) else "OPEN",
            "D2_disc_dressed_to_cutsum_same_scheme_comparator": "OPEN",
            "D3_uncertainty_propagation_from_dressed_backend": "OPEN",
        },
        "recommended_next_honest_step": {
            "id": "P2105_candidate",
            "goal": "build D2 comparator contract Disc_dressed vs CutSum on common basis with explicit row mapping",
        },
        "c3_gate_update": {
            "C3_d1_dressed_pole_residue_source_object": "COMPUTED" if (pre_ready and has_rows) else "OPEN",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": pre_ready,
            "d1_source_object_present": has_rows,
            "seed_level_only": True,
            "full_dressed_backend_import_executed": False,
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": False,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2104 S1054: strict D1 dressed-pole/residue source object",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Seed rows present: `{has_rows}`",
            f"- Seed row count: `{payload['d1_source_object']['rows_count']}`",
            "",
            "This stage exports a seed-level D1 source object for dressed pole/residue backend import.",
            "No full dressed Cutkosky closure claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
