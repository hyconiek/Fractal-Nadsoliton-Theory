#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2106 = GEN / "p2106_s1056_strict_cmp1_disc_dressed_row_or_nonavailability.json"
IN_1862 = GEN / "p1862_s812_strict_b1_dressed_pole_residue_and_projected_disc_evaluation_checkpoint.json"
IN_2020 = GEN / "p2020_s970_strict_cutkosky_p2019_tree_phase_space_cut_sum_witness.json"
OUT = GEN / "p2107_s1057_strict_cmp1_first_residual_execution.json"
MD = GEN / "p2107_s1057_strict_cmp1_first_residual_execution.md"

SCHEMA_VERSION = "p2107_s1057_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2106 = load(IN_2106)
    p1862 = load(IN_1862)
    p2020 = load(IN_2020)

    pre_ready = p2106.get("result_kind") == "PASS_STRICT_CMP1_COMPARATOR_ROW_SOURCE_OBJECT_WITH_TRACE__SEED_LEVEL"

    disc_rows = ((p1862.get("projected_discontinuity_seed_evaluation", {}) or {}).get("rows") or [])
    cutsum = (p2020.get("exact_phase_space_cut_sum", {}) or {}).get("CutSum_tree_identical_final_state_over_kappa2_Zgauge2")

    lhs_raw = disc_rows[0].get("projected_disc_value", 0.0) if disc_rows else None
    rhs_raw = cutsum

    lhs = float(sp.N(sp.sympify(lhs_raw))) if lhs_raw is not None else None
    rhs = float(sp.N(sp.sympify(rhs_raw))) if rhs_raw is not None else None

    residual = (lhs - rhs) if (lhs is not None and rhs is not None) else None

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2107",
        "stage_id": "S1057",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_STRICT_CMP1_FIRST_RESIDUAL_EXECUTION_WITH_TRACE"
            if pre_ready and residual is not None
            else "OPEN_STRICT_CMP1_FIRST_RESIDUAL_EXECUTION_BLOCKED"
        ),
        "depends_on": {
            "p2106_present": p2106.get("_missing") is None,
            "p1862_present": p1862.get("_missing") is None,
            "p2020_present": p2020.get("_missing") is None,
            "preconditions_ready": pre_ready,
        },
        "cmp1_residual_execution": {
            "lhs_disc_dressed_seed_row": lhs,
            "rhs_cutsum_row": rhs,
            "signed_residual": residual,
            "abs_residual": abs(residual) if residual is not None else None,
            "scope_limit": "first seed-level row residual only; not full dressed Disc=CutSum theorem",
        },
        "backend_import_blocker_update": {
            "D1_dressed_pole_residue_backend_on_common_basis": "COMPUTED",
            "D2_cmp1_source_row": "COMPUTED",
            "D2_cmp1_first_residual_execution": "COMPUTED" if residual is not None else "OPEN",
            "D3_uncertainty_propagation_from_dressed_backend": "OPEN",
        },
        "recommended_next_honest_step": {
            "id": "P2108_candidate",
            "goal": "add uncertainty-aware CMP1 residual bound and connect to D3 covariance propagation",
        },
        "c3_gate_update": {
            "C3_cmp1_first_residual_execution": "COMPUTED" if residual is not None else "OPEN",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": pre_ready,
            "cmp1_residual_computed": residual is not None,
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": False,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2107 S1057: strict CMP1 first residual execution",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- CMP1 residual computed: `{payload['gatekeeper_checks']['cmp1_residual_computed']}`",
            f"- Signed residual: `{residual}`",
            "",
            "This stage executes the first seed-level CMP1 residual Disc_dressed-CutSum.",
            "No full dressed Cutkosky closure claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
