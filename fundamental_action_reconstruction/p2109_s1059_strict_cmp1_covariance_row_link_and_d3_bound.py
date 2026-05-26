#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2108 = GEN / "p2108_s1058_strict_cmp1_uncertainty_aware_residual_bound.json"
IN_2107 = GEN / "p2107_s1057_strict_cmp1_first_residual_execution.json"
IN_2009 = GEN / "p2009_s959_strict_cutkosky_channelwise_tensor_coupled_covariance_classifier_witness.json"
OUT = GEN / "p2109_s1059_strict_cmp1_covariance_row_link_and_d3_bound.json"
MD = GEN / "p2109_s1059_strict_cmp1_covariance_row_link_and_d3_bound.md"

SCHEMA_VERSION = "p2109_s1059_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2108 = load(IN_2108)
    p2107 = load(IN_2107)
    p2009 = load(IN_2009)

    pre_ready = (
        p2108.get("result_kind") == "PASS_STRICT_CMP1_UNCERTAINTY_AWARE_BOUND_WITH_TRACE__SEED_LEVEL"
        and p2107.get("result_kind") == "PASS_STRICT_CMP1_FIRST_RESIDUAL_EXECUTION_WITH_TRACE"
    )

    residual = (p2107.get("cmp1_residual_execution", {}) or {}).get("signed_residual")
    residual_f = float(sp.N(sp.sympify(residual))) if residual is not None else None

    sigma_disc = float((p2108.get("cmp1_uncertainty_aware_residual_bound", {}) or {}).get("sigma_disc_seed_proxy", 0.0))

    C4 = np.array((p2009.get("backend_tensor_objects", {}) or {}).get("C4_channel_covariance", []), dtype=float)
    c4_valid = C4.shape == (4, 4)

    # Object-specific link for CMP1: first comparator row mapped to gg channel amplitude row.
    # This is a strict operational D3 row-link export, not a theorem-grade global transport proof.
    row_map = np.array([1.0, 0.0, 0.0, 0.0], dtype=float)
    sigma_cutsum = float(np.sqrt(max(0.0, row_map @ C4 @ row_map))) if c4_valid else 0.0

    sigma_total = float(np.sqrt(sigma_disc**2 + sigma_cutsum**2))
    if residual_f is not None:
        z95 = 1.96
        lo95 = residual_f - z95 * sigma_total
        hi95 = residual_f + z95 * sigma_total
    else:
        z95 = 1.96
        lo95 = None
        hi95 = None

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2109",
        "stage_id": "S1059",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_STRICT_CMP1_OBJECT_SPECIFIC_COVARIANCE_ROW_LINK_WITH_TRACE__D3_NONZERO_BOUND"
            if pre_ready and residual_f is not None and c4_valid and sigma_cutsum > 0.0
            else "OPEN_STRICT_CMP1_OBJECT_SPECIFIC_COVARIANCE_ROW_LINK_BLOCKED"
        ),
        "depends_on": {
            "p2108_present": p2108.get("_missing") is None,
            "p2107_present": p2107.get("_missing") is None,
            "p2009_present": p2009.get("_missing") is None,
            "preconditions_ready": pre_ready,
        },
        "cmp1_covariance_row_link": {
            "source_covariance_object": "p2009::backend_tensor_objects.C4_channel_covariance",
            "source_shape": list(C4.shape) if C4.size else [0, 0],
            "row_id": "CMP1",
            "row_map_to_channels": row_map.tolist(),
            "mapped_channel_interpretation": "CMP1 -> gg projected channel row",
            "sigma_cutsum_from_covariance_row_link": sigma_cutsum,
            "row_link_exported": c4_valid,
            "scope_limit": "object-specific row link only; not global D3 transport theorem",
        },
        "cmp1_uncertainty_aware_residual_bound_d3": {
            "residual_center": residual_f,
            "sigma_disc_seed_proxy": sigma_disc,
            "sigma_cutsum_from_row_link": sigma_cutsum,
            "sigma_total": sigma_total,
            "z95": z95,
            "residual_interval_95": [lo95, hi95],
            "contains_zero": (lo95 <= 0.0 <= hi95) if lo95 is not None else None,
        },
        "backend_import_blocker_update": {
            "D1_dressed_pole_residue_backend_on_common_basis": "COMPUTED",
            "D2_cmp1_source_row": "COMPUTED",
            "D2_cmp1_first_residual_execution": "COMPUTED",
            "D3_uncertainty_propagation_from_dressed_backend": (
                "COMPUTED_OBJECT_SPECIFIC_ROW_LINK"
                if c4_valid and sigma_cutsum > 0.0
                else "OPEN"
            ),
        },
        "recommended_next_honest_step": {
            "id": "P2110_candidate",
            "goal": "extend covariance row-link from CMP1-only to contract-level CMP-table and reconcile row-map with explicit dressed backend channel decomposition",
        },
        "c3_gate_update": {
            "C3_cmp1_first_residual_execution": "COMPUTED",
            "C3_cmp1_uncertainty_bound_seed_level": "COMPUTED",
            "C3_cmp1_object_specific_d3_row_link": "COMPUTED" if c4_valid and sigma_cutsum > 0.0 else "OPEN",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": pre_ready,
            "covariance_row_link_exported": c4_valid,
            "sigma_cutsum_nonzero": sigma_cutsum > 0.0,
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
            "# P2109 S1059: strict CMP1 object-specific covariance row link + D3 bound",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- sigma_cutsum (nonzero from row link): `{sigma_cutsum}`",
            f"- 95% interval: `{payload['cmp1_uncertainty_aware_residual_bound_d3']['residual_interval_95']}`",
            "",
            "This stage exports an object-specific covariance row link for CMP1 and recomputes the residual bound with nonzero dressed-side uncertainty.",
            "No global D3 transport theorem, full Cutkosky closure, or ToE closure is claimed.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
