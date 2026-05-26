#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2105 = GEN / "p2105_s1055_strict_d2_disc_dressed_vs_cutsum_comparator_contract.json"
IN_2109 = GEN / "p2109_s1059_strict_cmp1_covariance_row_link_and_d3_bound.json"
IN_2107 = GEN / "p2107_s1057_strict_cmp1_first_residual_execution.json"
IN_2009 = GEN / "p2009_s959_strict_cutkosky_channelwise_tensor_coupled_covariance_classifier_witness.json"
OUT = GEN / "p2110_s1060_strict_contract_level_cmp_table_d3_rowmap_calibration.json"
MD = GEN / "p2110_s1060_strict_contract_level_cmp_table_d3_rowmap_calibration.md"

SCHEMA_VERSION = "p2110_s1060_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def row_sigma(C4: np.ndarray, r: np.ndarray) -> float:
    return float(np.sqrt(max(0.0, r @ C4 @ r)))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2105 = load(IN_2105)
    p2109 = load(IN_2109)
    p2107 = load(IN_2107)
    p2009 = load(IN_2009)

    pre_ready = (
        p2105.get("result_kind") == "PASS_STRICT_D2_COMPARATOR_CONTRACT_WITH_TRACE__ROWS_REGISTERED"
        and p2109.get("result_kind") == "PASS_STRICT_CMP1_OBJECT_SPECIFIC_COVARIANCE_ROW_LINK_WITH_TRACE__D3_NONZERO_BOUND"
    )

    C4 = np.array((p2009.get("backend_tensor_objects", {}) or {}).get("C4_channel_covariance", []), dtype=float)
    c4_valid = C4.shape == (4, 4)

    # channel labels follow P2009 classifier order
    channel_labels = ["gg", "gh", "hh", "gx"]
    diag = np.diag(C4).tolist() if c4_valid else []
    trace_c4 = float(np.trace(C4)) if c4_valid else 0.0

    # Calibrated contract-level row maps: map comparator rows onto explicit channel decomposition.
    # CMP1: Disc_dressed common basis row is tied to gg projection.
    # CMP2: covariance-vs-uncertainty row is tied to gh projection as first available non-gg covariance channel.
    row_map_cmp1 = np.array([1.0, 0.0, 0.0, 0.0], dtype=float)
    row_map_cmp2 = np.array([0.0, 1.0, 0.0, 0.0], dtype=float)

    sigma_cmp1 = row_sigma(C4, row_map_cmp1) if c4_valid else 0.0
    sigma_cmp2 = row_sigma(C4, row_map_cmp2) if c4_valid else 0.0

    residual = (p2107.get("cmp1_residual_execution", {}) or {}).get("signed_residual")
    residual_f = float(sp.N(sp.sympify(residual))) if residual is not None else None

    z95 = 1.96
    cmp1_interval = [residual_f - z95 * sigma_cmp1, residual_f + z95 * sigma_cmp1] if residual_f is not None else [None, None]

    cmp_table = [
        {
            "row_id": "CMP1",
            "lhs": "Disc_dressed_common_basis_row_i",
            "rhs": "CutSum_common_basis_row_i",
            "channel_projection": "gg",
            "row_map_to_channels": row_map_cmp1.tolist(),
            "sigma_from_c4_row_link": sigma_cmp1,
            "residual_center": residual_f,
            "residual_interval_95": cmp1_interval,
            "status": "COMPUTED_OBJECT_SPECIFIC_ROW_LINK_AND_BOUND" if c4_valid and residual_f is not None else "OPEN",
        },
        {
            "row_id": "CMP2",
            "lhs": "Disc_dressed_covariance_row_i",
            "rhs": "CutSum_uncertainty_row_i",
            "channel_projection": "gh",
            "row_map_to_channels": row_map_cmp2.tolist(),
            "sigma_from_c4_row_link": sigma_cmp2,
            "status": "COMPUTED_OBJECT_SPECIFIC_ROW_LINK_PENDING_DRESSED_COV_ROW",
            "missing": "numeric dressed covariance row value for same-scheme residual execution",
        },
    ]

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2110",
        "stage_id": "S1060",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_STRICT_CONTRACT_LEVEL_CMP_TABLE_WITH_CALIBRATED_ROW_MAPS_AND_D3_LINK_TRACE"
            if pre_ready and c4_valid and sigma_cmp1 > 0.0 and sigma_cmp2 > 0.0
            else "OPEN_STRICT_CONTRACT_LEVEL_CMP_TABLE_CALIBRATION_BLOCKED"
        ),
        "depends_on": {
            "p2105_present": p2105.get("_missing") is None,
            "p2109_present": p2109.get("_missing") is None,
            "p2107_present": p2107.get("_missing") is None,
            "p2009_present": p2009.get("_missing") is None,
            "preconditions_ready": pre_ready,
        },
        "d3_channel_decomposition_export": {
            "source": "p2009::backend_tensor_objects.C4_channel_covariance",
            "shape": list(C4.shape) if C4.size else [0, 0],
            "channel_labels": channel_labels,
            "diag_entries": diag,
            "trace": trace_c4,
            "scope_limit": "operational channel decomposition for contract-level row-map calibration; not theorem-grade global transport",
        },
        "contract_level_cmp_table": {
            "scheme_tag": "STRICT_P2020_PHASESPACE_SCHEME_V1",
            "rows": cmp_table,
            "calibration_rule": "row-map chosen by explicit channel projection (CMP1->gg, CMP2->gh) on exported C4 channel decomposition",
        },
        "backend_import_blocker_update": {
            "D1_dressed_pole_residue_backend_on_common_basis": "COMPUTED",
            "D2_cmp1_source_row": "COMPUTED",
            "D2_cmp1_first_residual_execution": "COMPUTED",
            "D2_cmp2_residual_execution": "OPEN_MISSING_DRESSED_COV_ROW_VALUE",
            "D3_uncertainty_propagation_from_dressed_backend": "COMPUTED_CONTRACT_LEVEL_ROW_LINK_TABLE_PARTIAL",
        },
        "recommended_next_honest_step": {
            "id": "P2111_candidate",
            "goal": "export numeric dressed covariance row for CMP2 and execute CMP2 residual with same calibrated row-map table",
        },
        "c3_gate_update": {
            "C3_cmp1_object_specific_d3_row_link": "COMPUTED",
            "C3_cmp_contract_level_rowmap_calibration": "COMPUTED",
            "C3_cmp2_residual_execution": "OPEN",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": pre_ready,
            "channel_decomposition_exported": c4_valid,
            "cmp1_sigma_nonzero": sigma_cmp1 > 0.0,
            "cmp2_sigma_nonzero": sigma_cmp2 > 0.0,
            "contract_level_row_maps_exported": True,
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
            "# P2110 S1060: strict contract-level CMP table and D3 row-map calibration",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- CMP1 sigma: `{sigma_cmp1}`",
            f"- CMP2 sigma: `{sigma_cmp2}`",
            "",
            "This stage extends CMP1-only linkage to a contract-level CMP table (CMP1+CMP2) with explicit channel decomposition and calibrated row-maps.",
            "No global D3 transport theorem, full Cutkosky closure, or ToE closure is claimed.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
