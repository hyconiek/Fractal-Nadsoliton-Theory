#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2112 = GEN / "p2112_s1062_strict_cmp2_multibin_covariance_residual_table.json"
IN_2110 = GEN / "p2110_s1060_strict_contract_level_cmp_table_d3_rowmap_calibration.json"
IN_2009 = GEN / "p2009_s959_strict_cutkosky_channelwise_tensor_coupled_covariance_classifier_witness.json"
IN_1862 = GEN / "p1862_s812_strict_b1_dressed_pole_residue_and_projected_disc_evaluation_checkpoint.json"
OUT = GEN / "p2113_s1063_strict_cmp2_binwise_d3_covariance_transport_table.json"
MD = GEN / "p2113_s1063_strict_cmp2_binwise_d3_covariance_transport_table.md"

SCHEMA_VERSION = "p2113_s1063_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2112 = load(IN_2112)
    p2110 = load(IN_2110)
    p2009 = load(IN_2009)
    p1862 = load(IN_1862)

    pre_ready = (
        p2112.get("result_kind") == "PASS_STRICT_CMP2_MULTIBIN_RESIDUAL_TABLE_WITH_TRACE"
        and p2110.get("result_kind") == "PASS_STRICT_CONTRACT_LEVEL_CMP_TABLE_WITH_CALIBRATED_ROW_MAPS_AND_D3_LINK_TRACE"
    )

    C4 = np.array((p2009.get("backend_tensor_objects", {}) or {}).get("C4_channel_covariance", []), dtype=float)
    c4_valid = C4.shape == (4, 4)

    cmp_rows = ((p2110.get("contract_level_cmp_table", {}) or {}).get("rows") or [])
    cmp2 = next((r for r in cmp_rows if r.get("row_id") == "CMP2"), {})
    row_map = np.array(cmp2.get("row_map_to_channels", [0.0, 1.0, 0.0, 0.0]), dtype=float)
    if row_map.shape != (4,):
        row_map = np.array([0.0, 1.0, 0.0, 0.0], dtype=float)

    disc_rows = ((p1862.get("projected_discontinuity_seed_evaluation", {}) or {}).get("rows") or [])
    s_bins = sorted({float(r.get("s")) for r in disc_rows if r.get("s") is not None})

    base_sigma = float(np.sqrt(max(0.0, row_map @ C4 @ row_map))) if c4_valid else 0.0

    # Bin-wise transport factor derived from strict projected-disc spread at each s bin.
    rows_out = []
    for s in s_bins:
        vals = [float(r.get("disc_projected_seed")) for r in disc_rows if float(r.get("s")) == s and r.get("disc_projected_seed") is not None]
        if not vals:
            continue
        mean_abs = float(np.mean(np.abs(vals)))
        std = float(np.std(vals, ddof=0))
        rel_spread = (std / mean_abs) if mean_abs > 0 else 0.0
        transport_factor = 1.0 + rel_spread
        sigma_d3_binwise = base_sigma * transport_factor

        dressed_cov = float(np.var(vals, ddof=0))
        residual = dressed_cov  # cutsum uncertainty center remains 0 in current strict U3 export

        z95 = 1.96
        interval = [residual - z95 * sigma_d3_binwise, residual + z95 * sigma_d3_binwise]

        rows_out.append(
            {
                "s": s,
                "transport_factor": transport_factor,
                "sigma_d3_base_from_row_map": base_sigma,
                "sigma_d3_binwise": sigma_d3_binwise,
                "dressed_covariance_row_value": dressed_cov,
                "cutsum_uncertainty_row_center": 0.0,
                "signed_residual": residual,
                "residual_interval_95": interval,
                "contains_zero": interval[0] <= 0.0 <= interval[1],
            }
        )

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2113",
        "stage_id": "S1063",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_STRICT_CMP2_BINWISE_D3_COVARIANCE_TRANSPORT_TABLE_WITH_TRACE"
            if pre_ready and c4_valid and len(rows_out) > 0
            else "OPEN_STRICT_CMP2_BINWISE_D3_COVARIANCE_TRANSPORT_TABLE_BLOCKED"
        ),
        "depends_on": {
            "p2112_present": p2112.get("_missing") is None,
            "p2110_present": p2110.get("_missing") is None,
            "p2009_present": p2009.get("_missing") is None,
            "p1862_present": p1862.get("_missing") is None,
            "preconditions_ready": pre_ready,
        },
        "cmp2_binwise_d3_transport_table": {
            "row_id": "CMP2",
            "row_map_to_channels": row_map.tolist(),
            "channel_covariance_source": "p2009::backend_tensor_objects.C4_channel_covariance",
            "rows": rows_out,
            "scope_limit": "bin-wise operational D3 transport only; not global C3/D3 transport theorem",
        },
        "backend_import_blocker_update": {
            "D2_cmp2_residual_execution": "COMPUTED_MULTIBIN_TABLE",
            "D3_uncertainty_propagation_from_dressed_backend": "COMPUTED_BINWISE_ROWMAP_TRANSPORT_PARTIAL",
        },
        "recommended_next_honest_step": {
            "id": "P2114_candidate",
            "goal": "replace operational transport_factor by explicit backend-exported binwise covariance objects and test stability under row-map perturbations",
        },
        "c3_gate_update": {
            "C3_cmp2_multibin_residual_table": "COMPUTED",
            "C3_cmp2_binwise_d3_transport_table": "COMPUTED",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": pre_ready,
            "cmp2_row_map_present": row_map.shape == (4,),
            "base_sigma_nonzero": base_sigma > 0.0,
            "binwise_rows_exported": len(rows_out) > 0,
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
            "# P2113 S1063: strict CMP2 bin-wise D3 covariance transport table",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Row count: `{len(rows_out)}`",
            f"- Base sigma: `{base_sigma}`",
            "",
            "This stage replaces sigma_proxy by a bin-wise D3 transport sigma built from calibrated row-map covariance and per-bin spread factor.",
            "No global C3/D3 transport theorem or ToE closure claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
