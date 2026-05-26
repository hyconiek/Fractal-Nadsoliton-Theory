#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2110 = GEN / "p2110_s1060_strict_contract_level_cmp_table_d3_rowmap_calibration.json"
IN_1862 = GEN / "p1862_s812_strict_b1_dressed_pole_residue_and_projected_disc_evaluation_checkpoint.json"
IN_2101 = GEN / "p2101_s1051_strict_u3_residue_positivity_uncertainty_witness.json"
OUT = GEN / "p2111_s1061_strict_cmp2_covariance_row_export_and_residual_execution.json"
MD = GEN / "p2111_s1061_strict_cmp2_covariance_row_export_and_residual_execution.md"

SCHEMA_VERSION = "p2111_s1061_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2110 = load(IN_2110)
    p1862 = load(IN_1862)
    p2101 = load(IN_2101)

    pre_ready = p2110.get("result_kind") == "PASS_STRICT_CONTRACT_LEVEL_CMP_TABLE_WITH_CALIBRATED_ROW_MAPS_AND_D3_LINK_TRACE"

    disc_rows = ((p1862.get("projected_discontinuity_seed_evaluation", {}) or {}).get("rows") or [])
    first_s = disc_rows[0]["s"] if disc_rows else None
    same_s_values = [float(r["disc_projected_seed"]) for r in disc_rows if first_s is not None and r.get("s") == first_s]

    dressed_cov_row = float(np.var(same_s_values, ddof=0)) if same_s_values else None

    u3_rows = ((p2101.get("u3_uncertainty_table", {}) or {}).get("rows") or [])
    cutsum_unc_row = None
    if first_s is not None:
        for r in u3_rows:
            if r.get("s") == first_s:
                cutsum_unc_row = float(sp.N(sp.sympify(r.get("disc_uncertainty_abs", 0.0))))
                break

    cmp2_residual = (dressed_cov_row - cutsum_unc_row) if (dressed_cov_row is not None and cutsum_unc_row is not None) else None

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2111",
        "stage_id": "S1061",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_STRICT_CMP2_NUMERIC_COVARIANCE_ROW_EXPORT_AND_FIRST_RESIDUAL_WITH_TRACE"
            if pre_ready and cmp2_residual is not None
            else "OPEN_STRICT_CMP2_NUMERIC_COVARIANCE_ROW_EXPORT_BLOCKED"
        ),
        "depends_on": {
            "p2110_present": p2110.get("_missing") is None,
            "p1862_present": p1862.get("_missing") is None,
            "p2101_present": p2101.get("_missing") is None,
            "preconditions_ready": pre_ready,
        },
        "cmp2_covariance_row_execution": {
            "row_id": "CMP2",
            "dressed_covariance_row_value": dressed_cov_row,
            "dressed_covariance_row_source": "variance over p1862 projected_disc rows at fixed s (first available s)",
            "cutsum_uncertainty_row_value": cutsum_unc_row,
            "cutsum_uncertainty_row_source": "p2101 u3_uncertainty_table first matching s disc_uncertainty_abs",
            "signed_residual": cmp2_residual,
            "abs_residual": abs(cmp2_residual) if cmp2_residual is not None else None,
            "scope_limit": "first operational CMP2 row only; not full covariance transport theorem",
        },
        "backend_import_blocker_update": {
            "D1_dressed_pole_residue_backend_on_common_basis": "COMPUTED",
            "D2_cmp1_first_residual_execution": "COMPUTED",
            "D2_cmp2_residual_execution": "COMPUTED_FIRST_ROW",
            "D3_uncertainty_propagation_from_dressed_backend": "COMPUTED_CONTRACT_LEVEL_ROW_LINK_TABLE_PARTIAL",
        },
        "recommended_next_honest_step": {
            "id": "P2112_candidate",
            "goal": "execute multi-row CMP2 covariance residual table across all available s bins and attach explicit rowwise uncertainty intervals",
        },
        "c3_gate_update": {
            "C3_cmp_contract_level_rowmap_calibration": "COMPUTED",
            "C3_cmp2_first_residual_execution": "COMPUTED" if cmp2_residual is not None else "OPEN",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": pre_ready,
            "cmp2_dressed_cov_row_exported": dressed_cov_row is not None,
            "cmp2_cutsum_unc_row_linked": cutsum_unc_row is not None,
            "cmp2_residual_computed": cmp2_residual is not None,
            "full_d3_covariance_transport_proven": False,
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": False,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2111 S1061: strict CMP2 covariance row export and first residual execution",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- CMP2 dressed covariance row: `{dressed_cov_row}`",
            f"- CMP2 residual: `{cmp2_residual}`",
            "",
            "This stage exports the first numeric CMP2 dressed covariance row and executes the first CMP2 residual on the calibrated contract-level route.",
            "No global covariance transport theorem, full Cutkosky closure, or ToE closure is claimed.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
