#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2105 = GEN / "p2105_s1055_strict_d2_disc_dressed_vs_cutsum_comparator_contract.json"
IN_1862 = GEN / "p1862_s812_strict_b1_dressed_pole_residue_and_projected_disc_evaluation_checkpoint.json"
IN_2020 = GEN / "p2020_s970_strict_cutkosky_p2019_tree_phase_space_cut_sum_witness.json"
OUT = GEN / "p2106_s1056_strict_cmp1_disc_dressed_row_or_nonavailability.json"
MD = GEN / "p2106_s1056_strict_cmp1_disc_dressed_row_or_nonavailability.md"

SCHEMA_VERSION = "p2106_s1056_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2105 = load(IN_2105)
    p1862 = load(IN_1862)
    p2020 = load(IN_2020)

    pre_ready = p2105.get("result_kind") == "PASS_STRICT_D2_COMPARATOR_CONTRACT_WITH_TRACE__ROWS_REGISTERED"

    disc_rows = ((p1862.get("projected_discontinuity_seed_evaluation", {}) or {}).get("rows") or [])
    cutsum_obj = (p2020.get("exact_phase_space_cut_sum", {}) or {})

    cmp1_exported = bool(disc_rows) and bool(cutsum_obj)

    cmp1_object = {
        "row_id": "CMP1",
        "lhs_source": "p1862::projected_discontinuity_seed_evaluation.rows[0]" if disc_rows else "MISSING",
        "rhs_source": "p2020::exact_phase_space_cut_sum.CutSum_tree_identical_final_state_over_kappa2_Zgauge2" if cutsum_obj else "MISSING",
        "lhs_available": bool(disc_rows),
        "rhs_available": bool(cutsum_obj),
        "comparator_ready": cmp1_exported,
        "scope_limit": "seed-level comparator row only; not full dressed Disc=CutSum theorem",
    }

    result_kind = (
        "PASS_STRICT_CMP1_COMPARATOR_ROW_SOURCE_OBJECT_WITH_TRACE__SEED_LEVEL"
        if pre_ready and cmp1_exported
        else "OPEN_STRICT_CMP1_COMPARATOR_ROW_SOURCE_OBJECT_BLOCKED"
    )

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2106",
        "stage_id": "S1056",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": result_kind,
        "depends_on": {
            "p2105_present": p2105.get("_missing") is None,
            "p1862_present": p1862.get("_missing") is None,
            "p2020_present": p2020.get("_missing") is None,
            "preconditions_ready": pre_ready,
        },
        "cmp1_row_source_object": cmp1_object,
        "backend_import_blocker_update": {
            "D1_dressed_pole_residue_backend_on_common_basis": "COMPUTED",
            "D2_cmp1_source_row": "COMPUTED" if cmp1_exported else "OPEN",
            "D2_disc_dressed_to_cutsum_same_scheme_comparator": "OPEN_PENDING_EXECUTION",
            "D3_uncertainty_propagation_from_dressed_backend": "OPEN",
        },
        "recommended_next_honest_step": {
            "id": "P2107_candidate",
            "goal": "execute first numeric residual on CMP1 (Disc_dressed_row - CutSum_row) under same scheme or export object-specific nonavailability",
        },
        "c3_gate_update": {
            "C3_cmp1_row_source_object": "COMPUTED" if cmp1_exported else "OPEN",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": pre_ready,
            "cmp1_source_row_exported": cmp1_exported,
            "seed_level_only": True,
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": False,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2106 S1056: strict CMP1 Disc_dressed row source object",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- CMP1 source row exported: `{cmp1_exported}`",
            "",
            "This stage exports the first comparator row source object for D2 on seed level.",
            "No full dressed Cutkosky closure claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
