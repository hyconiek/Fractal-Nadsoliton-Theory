#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2104 = GEN / "p2104_s1054_strict_d1_dressed_pole_residue_source_object.json"
IN_2020 = GEN / "p2020_s970_strict_cutkosky_p2019_tree_phase_space_cut_sum_witness.json"
OUT = GEN / "p2105_s1055_strict_d2_disc_dressed_vs_cutsum_comparator_contract.json"
MD = GEN / "p2105_s1055_strict_d2_disc_dressed_vs_cutsum_comparator_contract.md"

SCHEMA_VERSION = "p2105_s1055_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2104 = load(IN_2104)
    p2020 = load(IN_2020)

    pre_ready = p2104.get("result_kind") == "PASS_STRICT_D1_DRESSED_POLE_RESIDUE_SOURCE_OBJECT_WITH_TRACE__SEED_LEVEL"

    comparator_rows = [
        {
            "row_id": "CMP1",
            "lhs": "Disc_dressed_common_basis_row_i",
            "rhs": "CutSum_common_basis_row_i",
            "status": "OPEN_MISSING_DRESSED_ROW",
            "missing": "dressed discontinuity row value on common basis",
        },
        {
            "row_id": "CMP2",
            "lhs": "Disc_dressed_covariance_row_i",
            "rhs": "CutSum_uncertainty_row_i",
            "status": "OPEN_MISSING_COVARIANCE_LINK",
            "missing": "same-scheme covariance propagation row",
        },
    ]

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2105",
        "stage_id": "S1055",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_STRICT_D2_COMPARATOR_CONTRACT_WITH_TRACE__ROWS_REGISTERED"
            if pre_ready
            else "OPEN_STRICT_D2_COMPARATOR_CONTRACT_BLOCKED"
        ),
        "depends_on": {
            "p2104_present": p2104.get("_missing") is None,
            "p2020_present": p2020.get("_missing") is None,
            "preconditions_ready": pre_ready,
        },
        "d2_comparator_contract": {
            "target": "Disc_dressed vs CutSum on common basis",
            "scheme_tag": "STRICT_P2020_PHASESPACE_SCHEME_V1",
            "rows": comparator_rows,
            "scope_limit": "contract registration only; comparator execution remains OPEN until dressed rows are exported",
        },
        "backend_import_blocker_update": {
            "D1_dressed_pole_residue_backend_on_common_basis": "COMPUTED",
            "D2_disc_dressed_to_cutsum_same_scheme_comparator": "OPEN_REGISTERED_CONTRACT",
            "D3_uncertainty_propagation_from_dressed_backend": "OPEN",
        },
        "recommended_next_honest_step": {
            "id": "P2106_candidate",
            "goal": "export first concrete Disc_dressed common-basis row (CMP1) or object-specific nonavailability for CMP1",
        },
        "c3_gate_update": {
            "C3_d2_comparator_contract": "COMPUTED",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": pre_ready,
            "comparator_rows_registered": len(comparator_rows) == 2,
            "d2_executed": False,
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": False,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2105 S1055: strict D2 Disc_dressed vs CutSum comparator contract",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Comparator rows registered: `{len(comparator_rows)}`",
            "",
            "This stage registers the same-scheme common-basis comparator contract rows for D2.",
            "No full dressed Cutkosky closure claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
